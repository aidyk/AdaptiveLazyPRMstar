# Adaptive Lazy PRMstar
A variation of Lazy PRM* to reduce the optimistic thrashing.

This work contains an implementaiton of adaptive lazy collision checking (edge) for Lazy PRMstar.
The main algorithmic framework follows Lazy PRMstar with DSPT, which was suggested by "Lazy Collision Checking in Asymptotically-Optimal Motion Planning" (Kris Hauser, 2015) and additionally applies a simple sample-based configuration-free space approximation for false negative error prediction.

This work is implemented with OMPL(Open Motion Planning Library, http://ompl.kavrakilab.org/), 
you can simply copy files into the "path-to-ompl/src/ompl/geometric/planners/prm/" and compile OMPL to run the code.

For more details, please visit our website (https://sglab.kaist.ac.kr/AdaptiveLazyCD/).
