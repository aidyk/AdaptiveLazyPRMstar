/*********************************************************************
  * Software License Agreement (BSD License)
  *
  *  Copyright (c) 2013, Willow Garage
  *  All rights reserved.
  *
  *  Redistribution and use in source and binary forms, with or without
  *  modification, are permitted provided that the following conditions
  *  are met:
  *
  *   * Redistributions of source code must retain the above copyright
  *     notice, this list of conditions and the following disclaimer.
  *   * Redistributions in binary form must reproduce the above
  *     copyright notice, this list of conditions and the following
  *     disclaimer in the documentation and/or other materials provided
  *     with the distribution.
  *   * Neither the name of Willow Garage nor the names of its
  *     contributors may be used to endorse or promote products derived
  *     from this software without specific prior written permission.
  *
  *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
  *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
  *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
  *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
  *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
  *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
  *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
  *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
  *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
  *  POSSIBILITY OF SUCH DAMAGE.
  *********************************************************************/

/* Author: A.I Dyk(Donghyuk Kim)*/

#include "ompl/geometric/planners/prm/AdaptivePRM.h"
#include "ompl/base/objectives/PathLengthOptimizationObjective.h"
#include "ompl/base/goals/GoalSampleableRegion.h"
#include "ompl/util/GeometricEquations.h"
#include "ompl/geometric/planners/prm/ConnectionStrategy.h"
#include "ompl/tools/config/SelfConfig.h"
#include "ompl/tools/debug/Profiler.h" // To enable and disable Profiling, please refer to the header file.
#include <boost/lambda/bind.hpp>
#include <boost/graph/astar_search.hpp>
#include <boost/graph/lookup_edge.hpp>
#include <boost/foreach.hpp>
#include <fstream>
#include <queue>
#include <utility>
#include <tuple>

#define foreach BOOST_FOREACH

#define BINARY
#define LESSLAZY
#define MORELAZY
#define LOGGING

namespace ompl {
namespace magic {
/** \brief The number of nearest neighbors to consider by
      default in the construction of the PRM roadmap */
static const unsigned int DEFAULT_NEAREST_NEIGHBORS_LAZY = 5;

/** \brief When optimizing solutions with lazy planners, this is the minimum
      number of path segments to add before attempting a new optimized solution
          extraction */
static const unsigned int MIN_ADDED_SEGMENTS_FOR_LAZY_OPTIMIZATION = 5;
}
}

ompl::geometric::AdaptivePRM::AdaptivePRM(const base::SpaceInformationPtr &si,
    bool starStrategy) :
  base::Planner(si, "AdaptivePRM"),
  starStrategy_(starStrategy),
  userSetConnectionStrategy_(false),
  maxDistance_(0.0),
  indexProperty_(boost::get(boost::vertex_index_t(), g_)),
  stateProperty_(boost::get(vertex_state_t(), g_)),
  radiusProperty_(boost::get(vertex_radius_t(), g_)),
  witnessProperty_(boost::get(vertex_witness_t(), g_)),
  costProperty_(boost::get(vertex_cost_t(), g_)),
  childrenProperty_(boost::get(vertex_children_t(), g_)),
  predecessorProperty_(boost::get(boost::vertex_predecessor_t(), g_)),
  colorProperty_(boost::get(boost::vertex_color_t(), g_)),
  weightProperty_(boost::get(boost::edge_weight_t(), g_)),
  vertexValidityProperty_(boost::get(vertex_flags_t(), g_)),
  edgeValidityProperty_(boost::get(edge_flags_t(), g_)),
  bestCost_(std::numeric_limits<double>::quiet_NaN()),
  iterations_(0),
  increaseIterations_(0),
  unPromisingCount_(0) {
  specs_.recognizedGoal = base::GOAL_SAMPLEABLE_REGION;
  specs_.approximateSolutions = false;
  specs_.optimizingPaths = true;

  Planner::declareParam<double>("range", this, &AdaptivePRM::setRange,
                                &AdaptivePRM::getRange, "0.:1.:10000.");

  if (!starStrategy_) {
    Planner::declareParam<unsigned int>("max_nearest_neighbors", this,
                                        &AdaptivePRM::setMaxNearestNeighbors, std::string("8:1000"));
  }

  addPlannerProgressProperty("iterations INTEGER",
                             boost::bind(&AdaptivePRM::getIterationCount, this));
  addPlannerProgressProperty("best cost REAL",
                             boost::bind(&AdaptivePRM::getBestCost, this));
  addPlannerProgressProperty("milestone count INTEGER",
                             boost::bind(&AdaptivePRM::getMilestoneCountString, this));
  addPlannerProgressProperty("edge count INTEGER",
                             boost::bind(&AdaptivePRM::getEdgeCountString, this));
}

ompl::geometric::AdaptivePRM::~AdaptivePRM() {
}

void ompl::geometric::AdaptivePRM::setup() {
  Planner::setup();
  tools::SelfConfig sc(si_, getName());
  sc.configurePlannerRange(maxDistance_);

  if (!nn_) {
    nn_.reset(tools::SelfConfig::getDefaultNearestNeighbors<Vertex>(this));
    nn_->setDistanceFunction(boost::bind(&AdaptivePRM::distanceFunction, this, _1, _2));
  }

  if (!connectionStrategy_) {
    if (starStrategy_) {
      connectionStrategy_ = KStarStrategy<Vertex>(boost::bind(
                              &AdaptivePRM::milestoneCount, this), nn_, si_->getStateDimension());
    } else {
      connectionStrategy_ = KBoundedStrategy<Vertex>
                            (magic::DEFAULT_NEAREST_NEIGHBORS_LAZY, maxDistance_, nn_);
    }
  }

  kPRMConstant_ = 1.1 * boost::math::constants::e<double>() + (boost::math::constants::e<double>() /
                  (double)si_->getStateDimension());
  double dimDbl = static_cast<double>(si_->getStateDimension());
  double prunedMeasure = si_->getSpaceMeasure();
  alpha_ = 0.8 * sqrt(static_cast<double>(si_->getStateDimension() * 0.5));

  // Setup optimization objective
  //
  // If no optimization objective was specified, then default to
  // optimizing path length as computed by the distance() function
  // in the state space.
  if (pdef_) {
    if (pdef_->hasOptimizationObjective()) {
      opt_ = pdef_->getOptimizationObjective();
    } else {
      opt_.reset(new base::PathLengthOptimizationObjective(si_));

      if (!starStrategy_) {
        opt_->setCostThreshold(opt_->infiniteCost());
      }
    }
  } else {
    OMPL_INFORM("%s: problem definition is not set, deferring setup completion...",
                getName().c_str());
    setup_ = false;
  }

  sampler_ = si_->allocStateSampler();

  // Get the measure of the entire space:
  prunedMeasure_ = si_->getSpaceMeasure();
  prunedCost_ = opt_->infiniteCost();
}

void ompl::geometric::AdaptivePRM::setRange(double distance) {
  maxDistance_ = distance;

  if (!userSetConnectionStrategy_) {
    connectionStrategy_.clear();
  }

  if (isSetup()) {
    setup();
  }
}

void ompl::geometric::AdaptivePRM::setMaxNearestNeighbors(unsigned int k) {
  if (starStrategy_) {
    throw Exception("Cannot set the maximum nearest neighbors for " + getName());
  }

  if (!nn_) {
    nn_.reset(tools::SelfConfig::getDefaultNearestNeighbors<Vertex>(this));
    nn_->setDistanceFunction(boost::bind(&AdaptivePRM::distanceFunction, this, _1, _2));
  }

  if (!userSetConnectionStrategy_) {
    connectionStrategy_.clear();
  }

  if (isSetup()) {
    setup();
  }
}

void ompl::geometric::AdaptivePRM::setProblemDefinition(
  const base::ProblemDefinitionPtr &pdef) {
  Planner::setProblemDefinition(pdef);
  clearQuery();
}

void ompl::geometric::AdaptivePRM::clearQuery() {
  startM_.clear();
  goalM_.clear();
  pis_.restart();
}

void ompl::geometric::AdaptivePRM::clear() {
  Planner::clear();
  freeMemory();

  if (nn_) {
    nn_->clear();
  }

  clearQuery();

  iterations_ = 0;
  bestCost_ = base::Cost(std::numeric_limits<double>::quiet_NaN());
  prunedCost_ = base::Cost(std::numeric_limits<double>::quiet_NaN());
  prunedMeasure_ = 0.0;
}

void ompl::geometric::AdaptivePRM::freeMemory() {
  foreach (Vertex v, boost::vertices(g_)) {
    si_->freeState(stateProperty_[v]);
  }

  g_.clear();
}

// Add newly sampled vertex and its adjancency edges connected to neigh neighbors.
ompl::geometric::AdaptivePRM::Vertex ompl::geometric::AdaptivePRM::addMilestone(
  base::State *state, bool isChecked) {
  Vertex m = boost::add_vertex(g_);
  stateProperty_[m] = state;
  radiusProperty_[m] = radiusProperty_[m] = std::numeric_limits<double>::infinity();
  witnessProperty_[m] = witnessProperty_[m] = nullptr;
  costProperty_[m] =  std::numeric_limits<double>::infinity();
  childrenProperty_[m] = new std::vector<Vertex>();
  predecessorProperty_[m] = nullptr;
  colorProperty_[m] = 0;
  vertexValidityProperty_[m] = (isChecked) ? VALIDITY_TRUE : VALIDITY_UNKNOWN;

  // Which milestones will we attempt to connect to?
  // const std::vector<Vertex> &neighbors = connectionStrategy_(m);

  std::vector<Vertex> neighbors;
  std::vector<double> neighbors_costs;
  const unsigned int max_number_of_neighbors = std::ceil(kPRMConstant_ * log(static_cast<double>(milestoneCount()) + 1u));

  neighbors.reserve(max_number_of_neighbors);
  neighbors_costs.reserve(max_number_of_neighbors);

  ompl::tools::Profiler::Begin("Nearest Neighbor Search");
  nn_->nearestK(m, max_number_of_neighbors, neighbors);
  ompl::tools::Profiler::End("Nearest Neighbor Search");

  foreach (Vertex n, neighbors) {
    if (witnessProperty_[n] != nullptr) {
      double dist = opt_->motionCost(witnessProperty_[n], stateProperty_[m]).value();
      // No directional information for radiusProperty...
      if (dist < radiusProperty_[m]) {
        witnessProperty_[m] = witnessProperty_[n];
        radiusProperty_[m] = dist;
      }
    }
  }

  foreach (Vertex n, neighbors) {
    // Symmetricity assumed here.
    const base::Cost weight = opt_->motionCost(stateProperty_[n], stateProperty_[m]);
    const Graph::edge_property_type properties(weight);

    // A.I : Adaptive Collision Checking
    if (!checkMotionAdaptively(n, m, weight.value())) {
      // Collision,
      continue;
    }

    // Collision-free,
    const Edge &e = boost::add_edge(n, m, properties, g_).first;
    edgeValidityProperty_[e] = VALIDITY_UNKNOWN;
  }
#ifdef BI_INHERIT
  if (witnessProperty_[m] != nullptr) {
    foreach (Vertex n, neighbors) {
      double dist = opt_->motionCost(witnessProperty_[m], stateProperty_[n]).value();

      if (dist < radiusProperty_[n]) {
        witnessProperty_[n] = witnessProperty_[m];
        radiusProperty_[n] = dist;
      }
    }
  }
#endif

  nn_->add(m);

  return m;
}

ompl::base::PlannerStatus ompl::geometric::AdaptivePRM::solve(const
    base::PlannerTerminationCondition &ptc) {
  // Initial checkup for start/goal configurations.
  checkValidity();

  base::GoalSampleableRegion *goal = dynamic_cast<base::GoalSampleableRegion *>
                                     (pdef_->getGoal().get());

  if (!goal) {
    OMPL_ERROR("%s: Unknown type of goal", getName().c_str());
    return base::PlannerStatus::UNRECOGNIZED_GOAL_TYPE;
  }

  // Add the valid start states as milestones
  while (const base::State *st = pis_.nextStart()) {
    Vertex st_vert = addMilestone(si_->cloneState(st));
    costProperty_[st_vert] = 0.0; // Initialize with 0 cost.
    startM_.push_back(st_vert);
  }

  if (startM_.size() == 0) {
    OMPL_ERROR("%s: There are no valid initial states!", getName().c_str());
    return base::PlannerStatus::INVALID_START;
  }

  if (!goal->couldSample()) {
    OMPL_ERROR("%s: Insufficient states in sampleable goal region",
               getName().c_str());
    return base::PlannerStatus::INVALID_GOAL;
  }

  // Ensure there is at least one valid goal state
  if (goal->maxSampleCount() > goalM_.size() || goalM_.empty()) {
    const base::State *st = goalM_.empty() ? pis_.nextGoal(ptc) : pis_.nextGoal();

    if (st) {
      goalM_.push_back(addMilestone(si_->cloneState(st)));
    }

    if (goalM_.empty()) {
      OMPL_ERROR("%s: Unable to find any valid goal states", getName().c_str());
      return base::PlannerStatus::INVALID_GOAL;
    }
  }

  unsigned long int nrStartStates = boost::num_vertices(g_);
  OMPL_INFORM("%s: Starting planning with %lu states already in datastructure",
              getName().c_str(), nrStartStates);

  ompl::tools::Profiler::Clear();
  ompl::tools::Profiler::Begin("Total");

  bestCost_ = opt_->infiniteCost();
  base::State *workState = si_->allocState();
  base::PathPtr bestSolution;
  bool fullyOptimized = false;
  bool someSolutionFound = false;

  /*
  int checkpoint = 1;
  auto t_start = std::clock();
  */

  // Grow roadmap in lazy fashion -- add edges without checking validity
  while (ptc == false) {
    ++iterations_;
    sampler_->sampleUniform(workState);
    /*
    double elapsed = (std::clock() - t_start) / CLOCKS_PER_SEC;
    if (elapsed >= checkpoint) {
      std::ofstream fileout;
      fileout.open("log2.txt", std::ofstream::out | std::ofstream::app);
      fileout << checkpoint << " " << bestCost_.value() << std::endl;
      fileout.close();
      checkpoint++;
    }*/

    // A.I : Pruning can be done here prior to checking its validity.
    /*
    if (!isPromising(workState)) {
      unPromisingCount_ += 1;
      continue;
    }*/

#ifdef MORELAZY
    // Add temporarily workState into Graph for NN, then remove again.
    Vertex m = boost::add_vertex(g_);
    stateProperty_[m] = workState;

    // Use nearest neighbor as a safty certificate in terms of collision.
    ompl::tools::Profiler::Begin("Nearest Neighbor Search");
    Vertex certificate = nn_->nearest(m), addedVertex;
    ompl::tools::Profiler::End("Nearest Neighbor Search");
    double safe_radius = radiusProperty_[certificate];
    double dist_to_certificate = distanceFunction(certificate, m); // Symmetry assumed.
    boost::clear_vertex(m, g_);
    boost::remove_vertex(m, g_); // Roll-back.

    if (dist_to_certificate < safe_radius) {
      // Add unknown vertex.
      addedVertex = addMilestone(si_->cloneState(workState), false);
    } else {
      ompl::tools::Profiler::Begin("Vertex Collision Checking");

      if (!si_->isValid(workState)) {
        ompl::tools::Profiler::End("Vertex Collision Checking");
        continue;
      }

      ompl::tools::Profiler::End("Vertex Collision Checking");

      // Add collision-free vertex.
      addedVertex = addMilestone(si_->cloneState(workState));
    }

#else
    // A.I : Yet wasting the collision information for sample right now...
    ompl::tools::Profiler::Begin("Vertex Collision Checking");

    if (!si_->isValid(workState)) {
      ompl::tools::Profiler::End("Vertex Collision Checking");
      continue;
    }

    ompl::tools::Profiler::End("Vertex Collision Checking");

    // Add collision-free vertices.
    Vertex addedVertex = addMilestone(si_->cloneState(workState), true);
#endif

    // DSPT update.
    Decrease(addedVertex);

    // Only support a single pair of start and goal node.
    Vertex startV = startM_[0];
    Vertex goalV = goalM_[0];
    base::PathPtr solution;

    do {
      if (predecessorProperty_[goalV] == nullptr) {
        break;
      }

      solution = constructSolution(startV, goalV);
    } while (!solution);

    if (solution) {
      someSolutionFound = true;
      base::Cost c = solution->cost(opt_);

      // A.I : If it is a non-optimal planner, it will be terminated at here,
      // will keep iterating otherwise.
      if (opt_->isSatisfied(c)) {
        fullyOptimized = true;
        bestSolution = solution;
        bestCost_ = c;
        break;
      } else {
        if (opt_->isCostBetterThan(c, bestCost_)) {
          bestSolution = solution;
          OMPL_INFORM("%.6lf -> %.6lf", bestCost_.value(), c.value());
          bestCost_ = c;
#ifdef PRUNETREE
          pruneTree(bestCost_);
#endif
        }
      }
    }
  }
  si_->freeState(workState);

  if (bestSolution) {
    base::PlannerSolution psol(bestSolution);
    psol.setPlannerName(getName());
    // If the solution was optimized, we mark it as such
    psol.setOptimized(opt_, bestCost_, fullyOptimized);
    pdef_->addSolutionPath(psol);
  }

  ompl::tools::Profiler::End("Total");

  /*
  std::ofstream fileout;
  fileout.open("log2.txt", std::ofstream::out | std::ofstream::app);
  fileout << checkpoint << " " << bestCost_.value() << std::endl << std::endl;
  fileout.close();
  */

  double avg_degree = 0.0;
  BGL_FORALL_VERTICES(vert, g_, Graph)
  avg_degree += boost::out_degree(vert, g_);
  avg_degree /= boost::num_vertices(g_);

  OMPL_INFORM("%s: Created %u vertices and %u edges with avg. degree %lf.", getName().c_str(),
              boost::num_vertices(g_) - nrStartStates, boost::num_edges(g_), avg_degree);

#ifdef LOGGING
  std::ofstream fileout;
  fileout.open("log.txt", std::ofstream::out | std::ofstream::app);
  fileout << bestCost_.value() << " " << boost::num_vertices(g_) <<
          " " << boost::num_edges(g_) << " " << avg_degree << std::endl;

  ompl::tools::Profiler::Status(fileout, true);
  fileout.close();

  OMPL_INFORM("%u", unPromisingCount_);
#endif

  return bestSolution ? base::PlannerStatus::EXACT_SOLUTION :
         base::PlannerStatus::TIMEOUT;
}

int ompl::geometric::AdaptivePRM::pruneTree(const base::Cost &pruneTreeCost) {
  return 0;
}

// A.I : Check edge collision checking from v1 to v2 and update wtness.
bool ompl::geometric::AdaptivePRM::checkMotion(const Vertex &v1, const Vertex &v2,
    const double weight) {
  ompl::tools::Profiler::ScopedBlock _profiler("Edge Collision Checking");
  /* Assume motion starts/ends in a valid configuration so v1/v2 are valid */
  bool result = true, witness_update = false;
  const base::State *s1 = stateProperty_[v1], *s2 = stateProperty_[v2];
  int nd = si_->getStateSpace()->validSegmentCount(s1, s2);
  double from_radius = (double)nd * radiusProperty_[v1] / weight;
  double to_radius = (double)nd * (1.0 - radiusProperty_[v2] / weight);

  if (nd > 1) {
    /* Temporary storage for the checked state */
    base::State *test = si_->allocState();
#ifndef BINARY
    // At j == 0 and nd -1 is corresponds to s1 and s2, respectively.
    for (int j = 1 ; j < nd ; j++) {
      si_->getStateSpace()->interpolate(s1, s2, (double)j / (double)nd, test);

      if (si_->isValid(test)) {
        continue;
      }
#ifdef LESSLAZY
      // A.I : Update witness and radius.
      // if ((double)j <= from_radius)  {
      if (radiusProperty_[v1] > ((double)j / (double)nd) * weight) {
        witness_update = true;
        witnessProperty_[v1] = test;
        radiusProperty_[v1] = ((double)j / (double)nd) * weight;
      }

      // if (to_radius <= (double)j) {
      if (radiusProperty_[v2] > (1.0 - (double)j / (double)nd) * weight) {
        witness_update = true;
        witnessProperty_[v2] = test;
        radiusProperty_[v2] = (1.0 - (double)j / (double)nd) * weight;
      }
#endif
      result = false;
      break;
    }
#else
    std::queue<std::pair<unsigned int, unsigned int> > q;
    q.push(std::make_pair(1, nd - 1);

    while (!q.empty()) {
      auto range = q.front();
      q.pop();
      unsigned int mid;

      mid = (range.first + range.second) / 2;
      si_->getStateSpace()->interpolate(s1, s2, (double)mid / (double)nd, test);

      if (!si_->isValid(test)) {
        result = false;
#ifdef LESSLAZY
        // A.I : Update witness and radius.
        double distance_to_obs = ((double)mid / (double)nd) * weight; // from v1
        if (radiusProperty_[v1] > distance_to_obs) {
          witness_update = true;
          witnessProperty_[v1] = test;
          radiusProperty_[v1] = distance_to_obs;
        }

        distance_to_obs = (1.0 - (double)mid / (double)nd) * weight; // from v2
        if (radiusProperty_[v2] > distance_to_obs) {
          witness_update = true;
          witnessProperty_[v2] = test;
          radiusProperty_[v2] = distance_to_obs;
        }
#endif
        break;
      }

      if (range.first < mid)
        q.push(std::make_pair(range.first, mid - 1));
      if (range.second < mid)
        q.push(std::make_pair(mid + 1, range.second));
      } // if mid == first, no more recursion.
    }

#endif

    if (!witness_update) {
      si_->freeState(test);  // If it is collision-free, no witness.
    }
  }

  // TODO : add counter.
  if (result) {
    // valid_++;
  } else {
    // invalid_++;
  }

  return result;
}

// A.I : Check adaptive edge collision checking from v1 to v2 and return witness closest from v1.
bool ompl::geometric::AdaptivePRM::checkMotionAdaptively(const Vertex &v1, const Vertex &v2,
                                                         const double weight) {
  ompl::tools::Profiler::ScopedBlock _profiler("Edge Collision Checking");
  /* Assume motion starts/ends in a valid configuration so v1/v2 are valid */
  const double r1 = radiusProperty_[v1], r2 = radiusProperty_[v2];
  bool result = true, witness_update = false;

  if (weight < r1 || weight < r2) {
    return true;
  }
  const base::State *s1 = stateProperty_[v1], *s2 = stateProperty_[v2];
  base::State *test = si_->allocState();

  // v1(from)-oriented
  si_->getStateSpace()->interpolate(s1, s2, r1 / weight, test); // x_proj
  double dist = opt_->motionCost(test, witnessProperty_[v1]).value();
  int isAligned = 0;
   if (dist < r1 * alpha_) { // IsAlign(.)
    isAligned = 1;
  } else {
    si_->getStateSpace()->interpolate(s2, s1, r2 / weight, test); // v2(to)-oriented
    dist = opt_->motionCost(test, witnessProperty_[v2]).value();
    if (dist < r2 * alpha_) { // IsAlign(.)
      isAligned = 2;
    }
  }

  if (isAligned) {
    double mid = r1 + (weight - r1 - r2) / 2.0;
    si_->getStateSpace()->interpolate(s1, s2, mid / weight, test); // x_mid
    // Test if they are aligned enough.
    if (!si_->isValid(test)) {
      result = false;
      if (mid < r1) { // Overlapping, both should update their witnesses.
        radiusProperty_[v1] = mid;
        radiusProperty_[v2] = weight - mid;
        witnessProperty_[v1] = witnessProperty_[v2] = test;
        witness_update = true;
      }
    }
  }

  if (!witness_update)
    si_->freeState(test);

  if (result) {
    // valid_++;
  } else {
    // invalid_++;
  }

  return result;
}

bool ompl::geometric::AdaptivePRM::isPromising(const base::State *s) {
  if (opt_->isFinite(bestCost_))
    return true;

  double dist = opt_->motionCost(stateProperty_[startM_[0]], s).value() +
                opt_->motionCost(s, stateProperty_[goalM_[0]]).value();
  return dist < bestCost_.value();
}

// outedge, inedge? - doesn't matter, need to scan all the neighbors.
void ompl::geometric::AdaptivePRM::Decrease(const Vertex &v) {
  ompl::tools::Profiler::ScopedBlock _profiler("Decrease");
  typedef std::pair<double, Vertex> weight_vertex;
  std::priority_queue<weight_vertex> pq;

  // Initialize cost of v, i.e., finding best parent vertex in G(g_).
  BGL_FORALL_OUTEDGES(v, e, g_, Graph) {
    Vertex w = target(e, g_);
    double weight = weightProperty_[e].value();

    if (costProperty_[v] > costProperty_[w] + weight) {
      predecessorProperty_[v] = w;
      costProperty_[v] = costProperty_[w] + weight;
    }
  }

  // No need to invoke cancelAdoption since v is newly sampled.
  if (predecessorProperty_[v] != nullptr) {
    childrenProperty_[predecessorProperty_[v]]->push_back(v);
  }

  // At this point, v has a best parent. From now on construct its subtree of descendants.

  pq.push(weight_vertex(-costProperty_[v], v)); // Invert the cost value for mimicking min-heap.

  while (!pq.empty()) {
    weight_vertex top = pq.top();
    pq.pop();

    double cost = -top.first; // Invert the cost value to be like min-heap.
    Vertex vert = top.second;

    if (cost > costProperty_[vert]) {
      continue;
    }

    BGL_FORALL_OUTEDGES(vert, e, g_, Graph) {
      Vertex w = target(e, g_);
      double weight = weightProperty_[e].value();
      double cost_w = costProperty_[w];

      if (cost_w > cost + weight) {
        costProperty_[w] = cost + weight;
        cancelAdoption(w);

        predecessorProperty_[w] = vert;
        childrenProperty_[vert]->push_back(w);
        pq.push(weight_vertex(-costProperty_[w], w));
      }
    }
  }

  // Now, DSPT is stable.
}

#define RED (increaseIterations_ + 1) // I know, it's bad #define.

void ompl::geometric::AdaptivePRM::Increase(const std::vector<Vertex> &vs) {
  ompl::tools::Profiler::ScopedBlock _profiler("Increase");
  // <Step 1. Preparation.
  // white is used for color of each vertex without initialization.
  // For each iteration, white is increased by 1, thus we can use it as
  // equal to or less than 'white' means 'white' color, 'red' otherwise.
  increaseIterations_ += 1;

  std::vector<Vertex> reds;
  typedef std::pair<double, Vertex> weight_vertex;
  std::priority_queue<weight_vertex> pq; //  Max-heap by default.

  for (unsigned int i = 0; i < vs.size(); i++) {
    pq.push(weight_vertex(-costProperty_[vs[i]], vs[i]));  // It works as if it is min-heap.
  }

  // <Step 2. Coloring
  while (!pq.empty()) {
    weight_vertex top = pq.top();
    pq.pop();

    double cost = -top.first;
    Vertex vert = top.second;

    if (cost > costProperty_[vert]) {
      continue;  // Instead of heap-improve
    }

    // If there exist a non-red neighbor q of z such that Dist(q) + w_(q, z) = D(z)
    // set pink, that means it can keep the current cost, thus it is not necessary to
    // iterate its children.
    // Otherwise, set red and enqueue all the children of z.
    bool pink_flag = false;
    BGL_FORALL_OUTEDGES(vert, e, g_, Graph) {
      Vertex w = target(e, g_);
      double weight = weightProperty_[e].value();

      if (colorProperty_[w] != RED && costProperty_[w] + weight == cost) {
        // Actually, '<' should not be happened all the time.
        // And even '==' would very rarely occur, but possible.
        pink_flag = true;

        cancelAdoption(vert);
        predecessorProperty_[vert] = w;
        childrenProperty_[w]->push_back(vert);
        break; // If there exsits, take anyone among them.
      }
    }

    if (pink_flag) {
      continue;
    }

    colorProperty_[vert] = RED; // Set to 'red'
    reds.push_back(vert);
    // Even with multiple starting red nodes, there will be no re-visit since each starting node is
    // a root node of sub'tree' in DSPT. That is, if statement within for loop might be useless.
    // Someone would be curious, e.g., then why do we have to use priority queue in step2 ?
    // Just for 'pink' case. Yeap. We need to identify all the other parent candidates are red or not
    // prior to checking current node.
    std::vector<Vertex> *children = childrenProperty_[vert];

    for (unsigned int i = 0; i < children->size(); i++) {
      pq.push(weight_vertex(-costProperty_[(*children)[i]], (*children)[i]));
    }
  }

  // 'pq' is empty at here.

  // <Step 3-a. Find best non-red parent for each red node.
  for (unsigned int i = 0; i < reds.size(); i++) {
    // TODO : need to be verified
    // Cost/predecessor initialization.
    costProperty_[reds[i]] = std::numeric_limits<double>::infinity();
    cancelAdoption(reds[i]);

    BGL_FORALL_OUTEDGES(reds[i], e, g_, Graph) {
      Vertex w = target(e, g_);
      double weight = weightProperty_[e].value();

      if (colorProperty_[w] == RED) {
        continue;  // If red, put aside for a while.
      }

      if (costProperty_[reds[i]] > costProperty_[w] + weight) {
        costProperty_[reds[i]] = costProperty_[w] + weight;
        predecessorProperty_[reds[i]] = w;
      }
    }

    if (predecessorProperty_[reds[i]] != nullptr) {
      childrenProperty_[predecessorProperty_[reds[i]]]->push_back(reds[i]);
    }

    pq.push(weight_vertex(-costProperty_[reds[i]], reds[i]));
  }

  // <Step 3-b. Propagate the changes; rewiring for 'red' nodes whether it can replace
  // existing parent node of near neighbors.
  while (!pq.empty()) {
    weight_vertex top = pq.top();
    pq.pop();

    double cost = -top.first;
    Vertex vert = top.second;

    if (costProperty_[vert] < cost) {
      continue;
    }

    BGL_FORALL_OUTEDGES(vert, e, g_, Graph) {
      Vertex w = target(e, g_);
      double weight = weightProperty_[e].value();

      if (colorProperty_[w] != RED) {
        continue;  // If not red, then skip.
      }

      if (cost + weight < costProperty_[w]) {
        costProperty_[w] = cost + weight;

        cancelAdoption(w);

        predecessorProperty_[w] = vert;
        childrenProperty_[vert]->push_back(w);
        pq.push(weight_vertex(-costProperty_[w], w));
      }
    }
  }

  // The end!
  // colorProperty_ doesn't need to be cleansed out since the increaing variable, 'RED'.
}

// TODO : sync between children & edges.
void ompl::geometric::AdaptivePRM::cancelAdoption(const Vertex &child) {
  ompl::tools::Profiler::ScopedBlock _profiler("Remove From Parent");

  if (predecessorProperty_[child] == nullptr) {
    return;
  }

  std::vector<Vertex> *children = childrenProperty_[predecessorProperty_[child]];

  if (children->empty()) {
    OMPL_DEBUG("ERROR");
  }

  for (unsigned int i = 0; i < children->size(); i++) if ((*children)[i] == child) {
    std::swap((*children)[i], children->back());
    children->pop_back();
    break;
  }

  predecessorProperty_[child] = nullptr;
}

// Vertex first.
ompl::base::PathPtr ompl::geometric::AdaptivePRM::constructSolution(
  const Vertex &start, const Vertex &goal) {
  ompl::tools::Profiler::ScopedBlock _profiler("Shortest Path Computation");
  std::vector<Vertex> solution_path;
  // ompl::RNG rng;

  // Construct a solution from DSPT.
  for (Vertex vert = goal; vert != nullptr; vert = predecessorProperty_[vert]) {
    solution_path.push_back(vert);
  }

  PathGeometric *p = new PathGeometric(si_);

  // Goal = 0, Start = n - 1.
  // TODO : From goal or start ? which one is better?

  auto from = solution_path.rbegin();

  if (from == solution_path.rend()) {
    return base::PathPtr();
  }

  for (auto to = from + 1; to != solution_path.rend(); ++to) {
    // std::vector<Vertex>::const_reverse_iterator from = solution_path.rbegin();
    // for (std::vector<Vertex>::const_reverse_iterator to = from + 1;
    // to != solution_path.rend(); ++to) {
    // OMPL_DEBUG(std::is_same<decltype(from), decltype(to)>::value ? "true" : "false");
#ifdef MORELAZY
    unsigned int &vd = vertexValidityProperty_[*to];
    const base::State *st = stateProperty_[*to];

    // Check collsion for unknown vertices.
    if ((vd & VALIDITY_TRUE) == 0) {
      ompl::tools::Profiler::Begin("Vertex Collision Checking");

      if (si_->isValid(st)) {
        vd |= VALIDITY_TRUE;
      }

      ompl::tools::Profiler::End("Vertex Collision Checking");
    }

    // If invalid...
    if ((vd & VALIDITY_TRUE) == 0) {
      // Nodes which lost the parent node i.e., 'to', becomes red!
      std::vector<Vertex> *children = childrenProperty_[*to], reds;

      for (unsigned int i = 0; i < children->size(); i++) {
        reds.push_back((*children)[i]);
        predecessorProperty_[(*children)[i]] = nullptr;
      }

      // Free from parent node, now it is in nowhere.
      cancelAdoption(*to);
      // Remove vertex from nearest neighbors data structure.
      nn_->remove(*to);
      // Free vertex state.
      si_->freeState(stateProperty_[*to]);
      // Remove all edges.
      boost::clear_vertex(*to, g_);
      // Remove the vertex.
      boost::remove_vertex(*to, g_);
      // Update DSPT
      Increase(reds);

      return base::PathPtr();
    }

#endif

    Edge e = boost::lookup_edge(*from, *to, g_).first; // Exhaustive search O(E) at worst case.
    unsigned int &evd = edgeValidityProperty_[e];

    if ((evd & VALIDITY_TRUE) == 0) { // Unknown
      bool result = true;
      result &= checkMotion(*from, *to, weightProperty_[e].value());

      if (result) {
        evd |= VALIDITY_TRUE;
      } else {
        boost::remove_edge(e, g_); // O(log(E/V)) time...
        cancelAdoption(*to);

        std::vector<Vertex> reds;
        reds.push_back(*to);
        Increase(reds);
        return base::PathPtr();
      }
    }

    from = to;
  }

  for (std::vector<Vertex>::const_reverse_iterator sol = solution_path.rbegin();
       sol != solution_path.rend(); ++sol) {
    p->append(stateProperty_[*sol]);
  }

  return base::PathPtr(p);
}

void ompl::geometric::AdaptivePRM::getPlannerData(base::PlannerData &data)
const {
  Planner::getPlannerData(data);

  // Explicitly add start and goal states. Tag all states known to be valid as 1.
  // Unchecked states are tagged as 0.
  for (size_t i = 0; i < startM_.size(); ++i) {
    data.addStartVertex(base::PlannerDataVertex(stateProperty_[startM_[i]], 1));
  }

  for (size_t i = 0; i < goalM_.size(); ++i) {
    data.addGoalVertex(base::PlannerDataVertex(stateProperty_[goalM_[i]], 1));
  }

  // Adding edges and all other vertices simultaneously
  foreach (const Vertex v, boost::vertices(g_)) {
    // const Vertex v1 = boost::source(e, g_);
    // const Vertex v2 = boost::target(e, g_);
    data.addEdge(base::PlannerDataVertex(stateProperty_[v]),
                 base::PlannerDataVertex(witnessProperty_[v]));

    // Add the reverse edge, since we're constructing an undirected roadmap
    // data.addEdge(base::PlannerDataVertex(stateProperty_[v2]),
    //             base::PlannerDataVertex(stateProperty_[v1]));

    float radius = static_cast<float>(radiusProperty_[v]);
    int disguised_float;
    memcpy(&disguised_float, &radius, sizeof(int));
    if (radius < 0.0f) printf("Invalid radius : %f\n", radius);

    // Add minimum radius to the closeset X_obs by exploiting tagState.
    data.tagState(stateProperty_[v], disguised_float);

    radius = -1.0f;
    memcpy(&disguised_float, &radius, sizeof(int));
    data.tagState(witnessProperty_[v], disguised_float);

    /*
    // Add tags for the newly added vertices
    data.tagState(stateProperty_[v1],
                  (vertexValidityProperty_[v1] & VALIDITY_TRUE) == 0 ? 0 : 1);
    data.tagState(stateProperty_[v2],
                  (vertexValidityProperty_[v2] & VALIDITY_TRUE) == 0 ? 0 : 1);
    */
  }
}
