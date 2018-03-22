# hello-world
This prog simulates the long rods problem with no intersection allowed in two dimensions on a square lattice. The rods are monodiepersed. It is a Monte Carlo(MC) simulation with heatbath dynamics where all the horizontal k-mers are evaporated and redeposited, then vertical k-mers and so on.  The initial conditions can be of two kinds: (1) all empty and (2) half vertical half horizontal. This is controlled using INITIAL_FLAG. The mean density, order parameter and second and fourth moment are measured. Horizontal rod is marked with 1 on lat while vertical rod is marked with 2.

The MC moves are nonlocal, thus it can equilibrate systems up to a very high density.
