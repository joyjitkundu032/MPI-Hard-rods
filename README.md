# Simulating a system on hard rods on a square lattice
This prog simulates a system of long rods with no intersection allowed in two dimensions on a square lattice in grand canonical ensemble. The chemical potential controls the density of the system. The rods are monodiepersed (of length k). It is a Monte Carlo(MC) simulation with heatbath dynamics where all the horizontal k-mers are evaporated and redeposited, then vertical k-mers and so on. The initial conditions can be of two kinds: (1) all empty and (2) half vertical half horizontal. This is controlled using INITIAL_FLAG. The mean density, order parameter, second and fourth moment along with the respective statistical error are measured. Horizontal rod is marked with 1 on lat while vertical rod is marked with 2.

The code requires two input files consisting of the list of probabilities for open and periodic boundary conditions. 

The MC moves are nonlocal, thus it can equilibrate systems up to a very high density that was inaccessible before. The code has been parallelized using MPI.
