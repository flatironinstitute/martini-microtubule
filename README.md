# martini-microtubule
Scripts for Iteratively Refined Distance Based Elastic Network adapted for microtubules.
1. iter0.py - Sets up an initial elastic network from relevant all-atom simulations. Uses distance-fluctuations to reduce the set of distance restraints, and scales the network strength on the basis of variance.
2. iterate,py - Iteratively optimizes the elastic network to reproduce target fluctuations.

## iter_0

### Inputs 
1. Coarse-grained frame of a protofilament.
2. Initial MARTINI molecule_0.itp file
3. Structure and Trajectory file from all-atom simulations.

### Outputs 
1. Distance fluctuation data from AA simulations in a pickled format.
2. Initial elastic-network files: enm and enm_a1b2

## Next Iterations

### Inputs 
1. Coarse-grained frame of a protofilament.
2. Initial MARTINI molecule_0.itp file
3. Structure and Trajectory file from previous coarse-grained simulations.

### Outputs 
1. elastic-network files for the new iteration: enm and enm_a1b2

