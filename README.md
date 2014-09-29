graphene
============

Implements grid matching of a template grid to an observed image


External dependencies
=============================
This package uses dm3lib_v099.py to read DM3 files. Download it from

https://code.google.com/p/diffraction-ring-profiler/source/browse/dm3lib_v099.py

and put it in the graphene/ folder.


Parameter tuning notes
=============================
In case you run the full pipeline and don't get the desired results, check the interemediate results before giving up! Also, use logging='DEBUG' in gridmatching.main() for troubleshooting.

Remember that no points are removed after the alternating graphcut procedure, they are only moved spatially.

Therefore check whether the local minima detection detects enough points; if it doesn't, focus on the parameters for that part. If too many points disappear due to the alternating graphcut, focus on the parameters for that part.

If the points are not positioned as you would expect, try tweaking the parameter 'beta' for the fine adjustment algorithm. In the extreme case of beta=0, there should be no spatial homogeneity - the points should end up in the closest local minimum. In the other extreme, they should position themselves such that the edge lengths are as homogeneous as possible. If you do not recognize this pattern from the experiments, make sure that a sufficient number of iterations are allowed for the simulated annealing, and perhaps play around with the temperature scheme.

Have fun!