/*********************************************************************
* Software License Agreement (BSD License)
*
*  Copyright (c) 2013, Rice University
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
*   * Neither the name of the Rice University nor the names of its
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

/* Author: A.I Dyk(Donghyuk Kim) */

#ifndef OMPL_GEOMETRIC_PLANNERS_PRM_ADAPTIVE_PRM_
#define OMPL_GEOMETRIC_PLANNERS_PRM_ADAPTIVE_PRM_

#include "ompl/geometric/planners/PlannerIncludes.h"
#include "ompl/datastructures/NearestNeighbors.h"
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/function.hpp>
#include <utility>
#include <vector>
#include <map>
#include <queue>

namespace ompl {
namespace base {
// Forward declare for use in implementation
OMPL_CLASS_FORWARD(OptimizationObjective);
}

namespace geometric {

/**
   @anchor gAdaptivePRM
   @par Short description
   AdaptivePRM is a planner that uses adaptive collision checking
   with approximate collision free space.
   @par External documentation
*/

/** \brief Lazy Probabilistic RoadMap planner */
class AdaptivePRM : public base::Planner {
 public:
  struct vertex_state_t {
    typedef boost::vertex_property_tag kind;
  };

  struct vertex_flags_t {
    typedef boost::vertex_property_tag kind;
  };

  struct vertex_radius_t {
    typedef boost::vertex_property_tag kind;
  };

  struct vertex_witness_t {
    typedef boost::vertex_property_tag kind;
  };

  struct vertex_cost_t {
    typedef boost::vertex_property_tag kind;
  };

  struct vertex_children_t {
    typedef boost::vertex_property_tag kind;
  };

  struct edge_flags_t {
    typedef boost::edge_property_tag kind;
  };

  /** @brief The type for a vertex in the roadmap. */
  typedef boost::adjacency_list_traits<boost::vecS, boost::listS,
          boost::undirectedS>::vertex_descriptor Vertex;

  /**
   @brief The underlying roadmap graph.

   @par Any BGL graph representation could be used here. Because we
   expect the roadmap to be sparse (m<n^2), an adjacency_list is more
   appropriate than an adjacency_matrix. We use listS for the vertex list
   because vertex descriptors are invalidated by remove operations if using vecS.

   @par Obviously, a ompl::base::State* vertex property is required.
   The incremental connected components algorithm requires
   vertex_predecessor_t and vertex_rank_t properties.
   If boost::vecS is not used for vertex storage, then there must also
   be a boost:vertex_index_t property manually added.

   @par Edges should be undirected and have a weight property.
   */
  typedef boost::adjacency_list <
  boost::vecS, boost::listS, boost::undirectedS,
        // Vertex properties.
        boost::property < vertex_state_t, base::State *,
        boost::property < boost::vertex_index_t, unsigned long int,
        boost::property < vertex_flags_t, unsigned int,
        boost::property < vertex_radius_t, double,
        boost::property < vertex_witness_t, base::State *,
        boost::property < vertex_cost_t, double,
        boost::property < vertex_children_t, std::vector<Vertex> *,
        boost::property < boost::vertex_color_t, unsigned int,
        boost::property < boost::vertex_predecessor_t, Vertex,
        boost::property < boost::vertex_rank_t, unsigned long int > > > > > > > > > >,
        // Edge properties.
        boost::property < boost::edge_weight_t, base::Cost,
        boost::property < edge_flags_t, unsigned int > >
        > Graph;

  /** @brief The type for an edge in the roadmap. */
  typedef boost::graph_traits<Graph>::edge_descriptor   Edge;

  /** @brief A nearest neighbors data structure for roadmap vertices. */
  typedef boost::shared_ptr< NearestNeighbors<Vertex> > RoadmapNeighbors;

  /** @brief A function returning the milestones that should be
   * attempted to connect to. */
  typedef boost::function<const std::vector<Vertex>&(const Vertex)> ConnectionStrategy;

  /** \brief Constructor */
  AdaptivePRM(const base::SpaceInformationPtr &si, bool starStrategy = false);

  virtual ~AdaptivePRM();

  /** \brief Set the maximum length of a motion to be added to the roadmap. */
  void setRange(double distance);

  /** \brief Get the range the planner is using */
  double getRange() const {
    return maxDistance_;
  }

  /** \brief Set a different nearest neighbors datastructure */
  template<template<typename T> class NN>
  void setNearestNeighbors() {
    nn_.reset(new NN<Vertex>());

    if (!userSetConnectionStrategy_) {
      connectionStrategy_.clear();
    }

    if (isSetup()) {
      setup();
    }
  }

  virtual void setProblemDefinition(const base::ProblemDefinitionPtr &pdef);

  /** \brief Set the connection strategy function that specifies the
   milestones that connection attempts will be make to for a
   given milestone.

   \par The behavior and performance of PRM can be changed drastically
   by varying the number and properties if the milestones that are
   connected to each other.

   \param pdef A function that takes a milestone as an argument and
   returns a collection of other milestones to which a connection
   attempt must be made. The default connection strategy is to connect
   a milestone's 10 closest neighbors.
   */
  void setConnectionStrategy(const ConnectionStrategy &connectionStrategy) {
    connectionStrategy_ = connectionStrategy;
    userSetConnectionStrategy_ = true;
  }

  /** \brief Convenience function that sets the connection strategy to the
   default one with k nearest neighbors.
   */
  void setMaxNearestNeighbors(unsigned int k);

  /** \brief Return the number of milestones currently in the graph */
  unsigned long int milestoneCount() const {
    return boost::num_vertices(g_);
  }

  /** \brief Return the number of edges currently in the graph */
  unsigned long int edgeCount() const {
    return boost::num_edges(g_);
  }

  virtual void getPlannerData(base::PlannerData &data) const;

  virtual void setup();

  virtual void clear();

  /** \brief Clear the query previously loaded from the ProblemDefinition.
      Subsequent calls to solve() will reuse the previously computed roadmap,
      but will clear the set of input states constructed by the previous call to solve().
      This enables multi-query functionality for AdaptivePRM. */
  void clearQuery();

  virtual base::PlannerStatus solve(const base::PlannerTerminationCondition &ptc);

  bool checkMotion(const Vertex &v1, const Vertex &v2, const double weight);
  bool checkMotionAdaptively(const Vertex &v1, const Vertex &v2, const double weight);

  void Decrease(const Vertex &v);
  void Increase(const std::vector<Vertex> &vs);

  bool isPromising(const base::State *s);
  int pruneTree(const base::Cost &pruneTreeCost);

  void cancelAdoption(const Vertex &v);

 protected:

  /** \brief Flag indicating validity of an edge of a vertex */
  static const unsigned int VALIDITY_UNKNOWN = 0;

  /** \brief Flag indicating validity of an edge of a vertex */
  static const unsigned int VALIDITY_TRUE    = 1;

  ///////////////////////////////////////
  // Planner progress property functions
  std::string getIterationCount() const {
    return boost::lexical_cast<std::string>(iterations_);
  }
  std::string getBestCost() const {
    return boost::lexical_cast<std::string>(bestCost_);
  }
  std::string getMilestoneCountString() const {
    return boost::lexical_cast<std::string>(milestoneCount());
  }
  std::string getEdgeCountString() const {
    return boost::lexical_cast<std::string>(edgeCount());
  }

  /** \brief Free all the memory allocated by the planner */
  void freeMemory();

  /** \brief Construct a milestone for a given state (\e state), store it in the nearest neighbors data structure
      and then connect it to the roadmap in accordance to the connection strategy. */
  Vertex addMilestone(base::State *state, bool isChecked = true);

  /** \brief Given two milestones from the same connected component, construct a path connecting them and set it as the solution */
  ompl::base::PathPtr constructSolution(const Vertex &start, const Vertex &goal);

  /** \brief Compute distance between two milestones (this is simply distance between the states of the milestones) */
  double distanceFunction(const Vertex a, const Vertex b) const {
    return si_->distance(stateProperty_[a], stateProperty_[b]);
  }

  /** \brief Given two vertices, returns a heuristic on the cost of the path connecting them.
      This method wraps OptimizationObjective::motionCostHeuristic */
  base::Cost costHeuristic(Vertex u, Vertex v) const;

  /** \brief Flag indicating whether the default connection strategy is the Star strategy */
  bool                                                   starStrategy_;

  /** \brief Function that returns the milestones to attempt connections with */
  ConnectionStrategy                                     connectionStrategy_;

  /** \brief Flag indicating whether the employed connection strategy was set by the user (or defaults are assumed) */
  bool
  userSetConnectionStrategy_;

  /** \brief The maximum length of a motion to be added to a tree */
  double                                                 maxDistance_;

  /** \brief Sampler user for generating random in the state space */
  base::StateSamplerPtr                                  sampler_;

  /** \brief Nearest neighbors data structure */
  RoadmapNeighbors                                       nn_;

  /** \brief Connectivity graph */
  Graph                                                  g_;

  /** \brief Array of start milestones */
  std::vector<Vertex>                                    startM_;

  /** \brief Array of goal milestones */
  std::vector<Vertex>                                    goalM_;

  /** \brief Access to the internal base::state at each Vertex */
  boost::property_map<Graph, boost::vertex_index_t>::type indexProperty_;

  /** \brief Access to the internal base::state at each Vertex */
  boost::property_map<Graph, vertex_state_t>::type       stateProperty_;

  /** \brief Access to the approximate collision-free radius at each Vertex */
  boost::property_map<Graph, vertex_radius_t>::type radiusProperty_;

  /** \brief Access to the approximate closest obstacle space at each Vertex */
  boost::property_map<Graph, vertex_witness_t>::type witnessProperty_;

  /** \brief Access to the cost from start configuration at each Vertex */
  boost::property_map<Graph, vertex_cost_t>::type costProperty_;

  /** \brief Access to the children of each Vertex */
  boost::property_map<Graph, vertex_children_t>::type childrenProperty_;

  /** \brief Access to the predecessor of each Vertex */
  boost::property_map<Graph, boost::vertex_predecessor_t>::type predecessorProperty_;

  /** \brief Access to the color of each Vertex used in Increase() */
  boost::property_map<Graph, boost::vertex_color_t>::type colorProperty_;

  /** \brief Access to the weights of each Edge */
  boost::property_map<Graph, boost::edge_weight_t>::type weightProperty_;

  /** \brief Access the validity state of a vertex */
  boost::property_map<Graph, vertex_flags_t>::type       vertexValidityProperty_;

  /** \brief Access the validity state of an edge */
  boost::property_map<Graph, edge_flags_t>::type         edgeValidityProperty_;

  /** \brief Objective cost function for PRM graph edges */
  base::OptimizationObjectivePtr                         opt_;

  base::Cost                                             bestCost_;

  unsigned long int                                      iterations_;

  /** \brief The number of 'Increase' invoked to avoid color initialization. */
  unsigned long int                                      increaseIterations_;

  double                                                 kPRMConstant_;
  double                                                 alpha_;

  /** brief A parameter for Collision Filter. */
  double                                                 kCFConstant_;
  base::Cost                                             prunedCost_;
  double                                                 prunedMeasure_;

  unsigned long int                                      unPromisingCount_;

  /** \brief Variables for custom benchmark. */
  unsigned int number_of_collision_checking_;
  unsigned int number_of_increase_;
};

}
}

#endif
