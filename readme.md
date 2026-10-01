# RayTraceHeatTransfer.jl

A Julia package for radiative heat transfer using quasi-Monte Carlo ray tracing and the Graph Equilibrium Radiative Transfer (GERT, see [Bielefeld, 2026](https://arxiv.org/abs/2512.22157)) methods. Solves grey, spectral and directional radiative equilibrium in 2D participating media and 3D surface enclosures, with unconditional energy conservation and exchange factor smoothing for machine precision reciprocity.

## Features

- **2D participating media** — absorbing, emitting, and scattering gases with enclosing surfaces
- **3D surface enclosures** — transparent media with semi-analytical view factors or ray tracing
- **Grey and spectral** — wavelength-independent or band-resolved radiation with automatic solver selection
- **Quasi-Monte Carlo sampling** — Sobol sequences by default, with results independent of the thread count and reproducible from a single seed
- **Adaptive spectral binning** — bins chosen automatically from κ(λ) samples to a user-set tolerance
- **Directional scattering and reflection** — Henyey–Greenstein or tabulated phase functions, specular or tabulated wall reflection, on angular bins with a user-set resolution (2D)
- **Unconditional energy conservation** — the GERT solve balances energy to machine precision for any number of rays, with or without smoothing
- **Exchange factor smoothing** — enforces reciprocity on the ray-traced factors to machine precision, which improves accuracy
- **Four-step workflow** — mesh → ray trace / view factors → smooth → solve
- **Plotting extensions** — GLMakie and Plots backends for mesh and field visualisation

## Installation

```julia
using Pkg
Pkg.add("RayTraceHeatTransfer")
```

For plotting, also install a Makie backend and/or Plots:

```julia
Pkg.add("GLMakie")   # for plotMesh (2D and 3D) and plotField (3D)
Pkg.add("Plots")     # for plotField (2D)
```

## Workflow

Every problem follows the same four steps:

1. **Mesh** — build the geometry and set boundary conditions.
2. **Exchange factors** — quasi-Monte Carlo ray tracing (2D participating media,
   3D surface enclosures) or analytical view factors (3D convex surface
   enclosures), producing exchange factor matrix `F_raw`.
3. **Smooth** — `smooth!(domain)` projects the exchange factor matrix `F_raw`
   onto the nearest matrix satisfying reciprocity and energy conservation,
   producing `F_smooth` and returning convergence diagnostics.
4. **Solve** — `solveEquilibrium!(domain, domain.F_smooth)` yields the
   steady state.

---

## Quickstart

A 1 × 1 m square filled with an absorbing gas (κ = 1 m⁻¹), the bottom wall at 1000 K and the other walls at 0 K, all black: mesh, trace, smooth and solve.

```julia
using RayTraceHeatTransfer
using GeometryBasics, StaticArrays

vertices   = SVector(Point2(0.0, 0.0), Point2(1.0, 0.0), Point2(1.0, 1.0), Point2(0.0, 1.0))   # a 1 m square
solidWalls = SVector(true, true, true, true)                       # all four walls are solid
face = PolyVolume2D{Float64}(vertices, solidWalls, 1, 1.0, 0.0)   # one spectral bin, κ = 1 m⁻¹, σₛ = 0
face.T_in_w  = [1000.0, 0.0, 0.0, 0.0]                             # bottom wall hot, the others cold [K]
face.epsilon = [1.0, 1.0, 1.0, 1.0]                                # black walls
face.T_in_g  = -1.0                                                # gas temperature unknown ...
face.q_in_g  = 0.0                                                 # ... in radiative equilibrium

mesh = RayTracingDomain2D([face], [(11, 11)])                      # 11 × 11 cells
mesh(10^7; method = :exchange)                                     # exchange factors by ray tracing
smooth!(mesh)                                                      # enforce reciprocity and energy conservation
solveEquilibrium!(mesh, mesh.F_smooth)                             # solve for the steady state

T_gas = [cell.T_g for cell in mesh.fine_mesh[1]]                   # gas temperature of every cell [K]
```

This is [Example 1](https://gert.net/examples/example-01.html) on gert.net, where the result is plotted and validated against Crosbie & Schrenker (1984).

## Examples

Worked examples, rendered with their code, figures and output, are on [gert.net](https://gert.net/examples/):

- [Example 1 — 2D Grey Participating Medium](https://gert.net/examples/example-01.html)
- [Example 2 — Spectral Greenhouse Atmosphere](https://gert.net/examples/example-02.html)
- [Example 3 — Line Spectrum vs Line-by-Line Reference](https://gert.net/examples/example-03.html)
- [Example 4 — Anisotropic Scattering vs a Discrete-Ordinates Reference](https://gert.net/examples/example-04.html)
- [Example 5 — Circular Enclosure from Triangular Elements](https://gert.net/examples/example-05.html)
- [Example 6 — 3D Surface Enclosure](https://gert.net/examples/example-06.html)
- [Example 7 — Triangulated Icosphere](https://gert.net/examples/example-07.html)
- [Example 8 — Mixed Triangular and Quadrilateral Faces](https://gert.net/examples/example-08.html)
- [Example 9 — Non-Convex Enclosure: an L-Shaped Duct](https://gert.net/examples/example-09.html)
- [Example 10 — Convergence in Rays and Mesh](https://gert.net/examples/example-10.html)
- [Example 11 — Large Scale: a Million Elements](https://gert.net/examples/example-11.html)

---

## References

The fundamental method is presented in:

- Bielefeld, N. M. (2026). A Radiation Exchange Factor Formulation with Proven Non-Negativity and Unconditional Energy Conservation. *arXiv preprint*, [arXiv:2512.22157](https://arxiv.org/abs/2512.22157).

It builds on the exchange factor method of:

- Larsen, M. E. (1983). *The Exchange Factor Method: An Alternative Zonal Formulation for Analysis of Radiating Enclosures Containing Participating Media*. PhD dissertation, The University of Texas at Austin. [OSTI](https://www.osti.gov/biblio/5004297)
- Larsen, M. E. & Howell, J. R. (1985). The exchange factor method: an alternative basis for zonal analysis of radiating enclosures. *Journal of Heat Transfer*, 107(4), 936–942. [doi:10.1115/1.3247524](https://doi.org/10.1115/1.3247524)

View factors and ray tracing:

- Narayanaswamy, A. (2015). An analytic expression for radiation view factor between two arbitrarily oriented planar polygons. *International Journal of Heat and Mass Transfer*, 91, 841–847. [doi:10.1016/j.ijheatmasstransfer.2015.07.131](https://doi.org/10.1016/j.ijheatmasstransfer.2015.07.131)
- Kerkhoff, J. A. & Wagner, M. J. (2021). A flexible thermal model for solar cavity receivers using analytical view factors. *Proceedings of the ASME 2021 15th International Conference on Energy Sustainability*, ES2021. [Link](https://asmedigitalcollection.asme.org/ES/proceedings-abstract/ES2021/84881/1114915)
- Kerkhoff, J. A. & Wagner, M. J. *viewFactor.m*. Energy Systems Optimization Lab, University of Wisconsin–Madison. [GitHub](https://github.com/uw-esolab/docs/tree/main/tools/viewfactor)
- Möller, T. & Trumbore, B. (1997). Fast, minimum storage ray-triangle intersection. *Journal of Graphics Tools*, 2(1), 21–28.
- Wald, I. (2007). On fast construction of SAH-based bounding volume hierarchies. *Proceedings of the 2007 IEEE Symposium on Interactive Ray Tracing*, 33–40.

Sampling:

- Sobol', I. M. (1967). On the distribution of points in a cube and the approximate evaluation of integrals. *USSR Computational Mathematics and Mathematical Physics*, 7(4), 86–112.
- Joe, S. & Kuo, F. Y. (2008). Constructing Sobol sequences with better two-dimensional projections. *SIAM Journal on Scientific Computing*, 30(5), 2635–2654.

Smoothing:

- Dykstra, R. L. (1983). An algorithm for restricted least squares regression. *Journal of the American Statistical Association*, 78(384), 837–842.
- Larsen, M. E. & Howell, J. R. (1986). Least-squares smoothing of direct-exchange areas in zonal analysis. *Journal of Heat Transfer*, 108(1), 239–242.

Physics and validation references used in the examples:

- Crosbie, A. L. & Schrenker, R. G. (1984). Radiative transfer in a two-dimensional rectangular medium exposed to diffuse radiation. *Journal of Quantitative Spectroscopy and Radiative Transfer*, 31(4), 339–372.
- Heaslet, M. A. & Warming, R. F. (1965). Radiative transport and wall temperature slip in an absorbing planar medium. *International Journal of Heat and Mass Transfer*, 8(7), 979–994. [doi:10.1016/0017-9310(65)90083-9](https://doi.org/10.1016/0017-9310(65)90083-9)
- Henyey, L. G. & Greenstein, J. L. (1941). Diffuse radiation in the galaxy. *The Astrophysical Journal*, 93, 70–83. [doi:10.1086/144246](https://doi.org/10.1086/144246)
- Howell, J. R., Mengüç, M. P., Daun, K. & Siegel, R. (2021). *Thermal Radiation Heat Transfer*, 7th ed. CRC Press.
- Modest, M. F. & Mazumder, S. (2022). *Radiative Heat Transfer*, 4th ed. Academic Press. [doi:10.1016/C2018-0-03206-5](https://doi.org/10.1016/C2018-0-03206-5)

Software:

- Bezanson, J., Edelman, A., Karpinski, S. & Shah, V. B. (2017). Julia: a fresh approach to numerical computing. *SIAM Review*, 59(1), 65–98.
- [Sobol.jl](https://github.com/JuliaMath/Sobol.jl): Sobol sequences in Julia.
- [ConvolutionInterpolations.jl](https://github.com/NikoBiele/ConvolutionInterpolations.jl): interpolation used for the spectral binning and the Planck tables.
- Danisch, S. & Krumbiegel, J. (2021). Makie.jl: flexible high-performance data visualization for Julia. *Journal of Open Source Software*, 6(65), 3349.
- Christ, S., Schwabeneder, D., Rackauckas, C., Borregaard, M. K. & Breloff, T. (2023). Plots.jl — a user extendable plotting API for the Julia programming language. *Journal of Open Research Software*, 11(1), 5.

## Authors

The primary author, developer and maintainer of this repository is Nikolaj Maack Bielefeld.

The functions for calculating 3D view factors analytically were originally written for MATLAB by Jacob A. Kerkhoff and Michael J. Wagner of University of Wisconsin-Madison, Energy Systems Optimization Lab, as described in [Kerkhoff & Wagner (2021)](https://asmedigitalcollection.asme.org/ES/proceedings-abstract/ES2021/84881/1114915).

## Declaration of AI Assistance

RayTraceHeatTransfer.jl has been developed since 2024 as a collaboration between the author and Claude, Anthropic's AI model, across its successive versions. The work has been shared throughout: ideas were proposed and challenged in both directions, algorithms were designed and implemented together, and much of the code, the tests, the validation references and the documentation took shape in that collaborative spirit. Most of the package's key ideas emerged in these discussions.

The author directed the project, made the final decisions, and reviewed and applied every change himself, and he takes full responsibility for the code, the methods and the scientific content. Everything is validated independently of how it was written: against analytical limits and published reference solutions, as shown in the examples on [gert.net](https://gert.net/examples/), and by the package test suite.