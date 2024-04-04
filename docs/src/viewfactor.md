### View Factors in 3D

As a utility functionality this repository also features a function to calculate view factors between arbitrary polygons separated by a fully transparent medium.
It is the long term intention to make RayTraceHeatTransfer.jl into a 3D ray tracing code and by then the 3D view factor function will be valuable for validation.

The original authors of the MATLAB-implementation of this function are: Jacob A. Kerkhoff and Michael J. Wagner of University of Wisconsin-Madison, Energy Systems Optimization Lab. The MATLAB implementation is available at the following repository: https://github.com/uw-esolab/docs/tree/main/tools/viewfactor .
This function was originally written for the paper: "A Flexible Thermal Model for Solar Cavity Receivers Using Analytical View Factors" https://doi.org/10.1115/ES2021-63810 .

The function was translated to the Julia programming language by the author of RayTraceHeatTransfer.jl during studies at Aalborg University, Thermal Energy and Process Engineering.
The ray tracing part of the original function was not translated. Published in RayTraceHeatTransfer.jl with permission from the original authors.

The analytical solution was derived in the following paper: Narayanaswamy, Arvind. "An analytic expression for radiation view factor between two arbitrarily oriented planar polygons." International Journal of Heat and Mass Transfer 91 (2015): 841-847 https://doi.org/10.1016/j.ijheatmasstransfer.2015.07.131 .

viewFactor(POLY_A, POLY_B) analytically computes the view factor from POLY_A to POLY_B and from POLY_B to POLY_A. Input arguments are in the form of 3x2 or Nx3 arrays, where each row corresponds to a vertex of the polygon, the columns refer to the X,Y,Z coordinates of the vertices. If Z coordinates are omitted, they are assumed to be zero.

The following must be true for both polygons:
1) polygons are planar (all vertices lie in the same plane)
2) polygons are simple (no self-intersecting polygons)
3) polygons are convex (in theory, concave polygons should work, but
    this remains untested)

Below an example is given:

```julia
POLY_A = [0.0 0.0 0.0
        1.0 0.0 0.0
        1.0 1.0 0.0
        0.0 1.0 0.0]
POLY_B = [0.0 0.0 1.0
        1.0 0.0 1.0
        1.0 1.0 1.0
        0.0 1.0 1.0]
resultPaper = 0.199825
F_AB, F_BA, area_A, area_B = viewFactor(POLY_A, POLY_B)
```