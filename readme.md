# RayTraceHeatTransfer.jl

A Julia package for radiative heat transfer using Monte Carlo ray tracing and the Graph Equilibrium Radiative Transfer (GERT) methods. Solves grey and spectral radiative equilibrium in 2D participating media and 3D surface enclosures, with exchange factor smoothing for machine precision energy-conserving solutions.

## Features

- **2D participating media** — absorbing, emitting, and scattering gases with enclosing surfaces
- **3D surface enclosures** — transparent media with analytical view factors (Narayanaswamy 2015)
- **Grey and spectral** — wavelength-independent or band-resolved radiation with automatic solver selection
- **Exchange factor smoothing** — iterative reciprocity enforcement for machine-precision energy conservation
- **Callable domains** — `mesh(N_rays; method=:exchange)` runs ray tracing directly on the domain object
- **Plotting extensions** — GLMakie and Plots backends for mesh and field visualisation

## Installation

```julia
using Pkg
Pkg.add("RayTraceHeatTransfer")
```

For plotting, also install a Makie backend and/or Plots:

```julia
Pkg.add("GLMakie")   # for plotMesh (2D and 3D)
Pkg.add("Plots")     # for plotField (2D)
```

---

## Example 1 — 2D Grey Participating Medium

This example solves radiative equilibrium in a 1 × 1 m square enclosure filled with an absorbing gas (κ = 1 m⁻¹, no scattering). The bottom wall is held at 1000 K and all other walls are at 0 K; all surfaces are black (ε = 1). The gas temperature field is found by solving the GERT system after computing exchange factors with Monte Carlo ray tracing.

### Step 1: Define the geometry and mesh

```julia
using RayTraceHeatTransfer
using GeometryBasics, StaticArrays

Ndim = 11  # 11 × 11 elements

vertices = SVector(
    Point2(0.0, 0.0),
    Point2(1.0, 0.0),
    Point2(1.0, 1.0),
    Point2(0.0, 1.0)
)
solidWalls = SVector(true, true, true, true)

face = PolyVolume2D{Float64}(vertices, solidWalls, 1, 1.0, 0.0);  # κ=1, σₛ=0

face.T_in_w  = [1000.0, 0.0, 0.0, 0.0]   # bottom hot, rest cold
face.epsilon = [1.0, 1.0, 1.0, 1.0]       # black walls
face.T_in_g  = -1.0                         # unknown (solve for this)
face.q_in_g  = 0.0                          # radiative equilibrium

mesh = RayTracingDomain2D([face], [(Ndim, Ndim)]);
```

### Understanding the mesh numbering

Each gas volume and solid wall surface in the mesh receives a global index that corresponds to a row/column in the exchange factor matrix. Use `plotMesh` with the `volumeNumbers` and `wallNumbers` keyword arguments to visualise specific indices:

```julia
using GLMakie

fig = Figure(size = (700, 700))
ax  = Axis(fig[1, 1], aspect = DataAspect(), xlabel = "x (m)", ylabel = "y (m)",
           title = "Mesh with element numbering (11 × 11)")

# Show volume indices along the vertical centerline
center_col = div(Ndim + 1, 2)
centerline_vols = [center_col + (row - 1) * Ndim for row in 1:Ndim]
bottom_wall_indices = [1; collect(3:Ndim+1)]

plotMesh(ax, mesh)
plotMesh(ax, mesh; volumeNumbers = centerline_vols)
plotMesh(ax, mesh; wallNumbers = bottom_wall_indices)

fig
save("fig/mesh_numbering.png", fig, px_per_unit=3)
```

![Mesh numbering](fig/mesh_numbering.png)

Volume elements are labelled **g*i*** and wall surfaces **w*i***. The indices shown here are the same indices used in the exchange factor matrix `mesh.F_smooth` and in the system matrices returned by `buildSystemMatrices!`.

### Step 2: Ray trace and solve

```julia
mesh(10_000_000; method = :exchange)          # Monte Carlo ray tracing
solveEquilibrium!(mesh, mesh.F_smooth)       # solve GERT system
```

### Step 3: Validate against Crosbie & Schrenker (1984)

The analytical solution for the dimensionless source function S(τ) = (T/T_hot)⁴ along the centerline of this problem is given by Crosbie & Schrenker (1984). Extracting the centerline temperatures from the solved mesh and comparing:

```julia
using Plots

# --- Left panel: solution temperature field via plotField ---
p1 = plotField(mesh; field = :T)

# Extract centerline temperatures
all_temps  = [ff.T_g for ff in mesh.fine_mesh[1]]
Tg_matrix  = reshape(all_temps, Ndim, Ndim)
centerline = Tg_matrix[div(Ndim + 1, 2), :]

# Dimensionless source function
S_computed = (centerline ./ 1000.0) .^ 4
tau_centers = range(1 / (2Ndim), 1 - 1 / (2Ndim), length = Ndim)

# --- Crosbie & Schrenker (1984) analytical reference ---
tau_ref = [0.0, 0.00611, 0.02037, 0.04251, 0.07216, 0.10884, 0.15194,
           0.20076, 0.25449, 0.31225, 0.37309, 0.43602, 0.50000, 0.56398,
           0.62691, 0.68775, 0.74551, 0.79924, 0.84806, 0.89116, 0.92784,
           0.95749, 0.97963, 0.99390, 1.00000]

S_ref  = [0.6293, 0.6198, 0.6017, 0.5767, 0.5460, 0.5108, 0.4724,
          0.4323, 0.3919, 0.3525, 0.3153, 0.2810, 0.2500, 0.2224,
          0.1981, 0.1768, 0.1584, 0.1424, 0.1287, 0.1171, 0.1073,
          0.0992, 0.0930, 0.0885, 0.0863]

# --- Right panel: centerline validation ---
p2 = Plots.plot(tau_ref, S_ref,
    linewidth = 2, color = :black, label = "Analytical (C & S)",
    xlabel = "Optical depth τ",
    ylabel = "Dimensionless source function S(τ)",
    title = "Centerline validation",
    legend = :topright,
    dpi = 1000,
    guidefontsize = 12,
    tickfontsize = 10)

Plots.scatter!(p2, tau_centers, S_computed,
    color = :dodgerblue, markersize = 5, label = "RayTraceHeatTransfer.jl")
Plots.plot!(p2, top_margin=8Plots.mm, bottom_margin=8Plots.mm)

# --- Combined figure ---
Plots.plot!(p1, guidefontsize=12, tickfontsize=10, 
            left_margin=5Plots.mm, right_margin=10Plots.mm)
p = Plots.plot(p1, p2, layout = (2,1), size = (600, 800), dpi=1000)
display(p)
Plots.savefig(p, "fig/validation_2d_grey.png")
```

![2D grey validation](fig/validation_2d_grey.png)

The top panel shows the 2D temperature field; the bottom panel compares the computed centerline source function (blue dots) with the analytical reference (black line). Agreement is excellent for 10⁷ rays on an 11 × 11 mesh.

### References

> Crosbie, A. L. & Schrenker, R. G. (1984). "Radiative transfer in a two-dimensional rectangular medium exposed to diffuse radiation." *Journal of Quantitative Spectroscopy and Radiative Transfer*, 31(4), 339–372.

> Bielefeld, N. M. (2025). "Graph Equilibrium Radiative Transfer." *arXiv preprint*, [arXiv:2512.22157](https://arxiv.org/abs/2512.22157).

---

## Example 2 — Spectral Greenhouse Atmosphere

This example models a simplified planetary atmosphere to demonstrate the spectral solver. The atmosphere is transparent in the visible and opaque in the infrared, producing a greenhouse effect: solar radiation penetrates to the surface, while thermal emission from the warm surface is trapped by the absorbing gas. The equilibrium surface temperature — not prescribed — emerges far above the bare blackbody value.

The geometry is a vertical stack of 20 sub-enclosures representing atmospheric layers, each with spectrally varying absorption. A thin volume at the top emits at solar temperature, acting as the radiation source. The domain is wide relative to its height, approximating a 1D atmosphere.

### Step 1: Define the atmosphere

```julia
using RayTraceHeatTransfer
using GeometryBasics, StaticArrays

n_bins       = 40            # spectral bins
atm_height   = 100_000.0     # atmosphere height (m)
L            = atm_height    # normalization length
N_layers     = 20            # atmospheric layers
width        = 100.0         # normalized width (wide domain ≈ 1D)
scale_height = 15_900.0      # density scale height (m)
T_sun        = 5800.0        # solar temperature (K)
q_solar      = 2 * 2600.0    # isotropic solar flux (W/m²)
κ_vis        = 0.01          # visible absorption coefficient
κ_ir         = 100.0         # infrared absorption coefficient
λ_min        = 1e-9          # minimum wavelength (m)
λ_max        = 1.0           # maximum wavelength (m)
stretch      = 5.0           # spatial layer clustering near surface
```

The spectral range spans from 1 nm to 1 m — wide enough to capture the full Planck distribution at all temperatures in the problem. An insufficient spectral range forces energy into edge bins and degrades the solution.

### Step 2: Build the spectral grid and layer geometry

```julia
# Spectral grid (log-spaced in wavelength)
λ_edges  = 10 .^ range(log10(λ_min), log10(λ_max), length = n_bins + 1)
λ_center = [sqrt(λ_edges[b] * λ_edges[b+1]) for b in 1:n_bins]

# Log-spaced spatial layers: thin near the surface where temperature
# gradients are steepest, thick higher up where the atmosphere thins
layer_param      = range(0.0, 1.0, length = N_layers + 1)
layer_edges_norm = [(exp(stretch * t) - 1) / (exp(stretch) - 1) for t in layer_param]

# Solar volume: a thin layer at the top whose emission matches the desired
# irradiance. This avoids modifying the solver for spectral boundary fluxes.
sun_layer_height = 1000.0     # 1 km thick
κ_sun = q_solar * L / (4 * 5.670374419e-8 * T_sun^4 * sun_layer_height)

normalized_scale_height = scale_height / L
```

### Step 3: Assemble the atmospheric layers

Each layer has a spectrally varying absorption coefficient: a sigmoid transition around λ = 4 μm separates the transparent visible window (κ ≈ 0.01) from the opaque infrared (κ ≈ 100), scaled by the local atmospheric density. This spectral asymmetry is the mechanism behind the greenhouse effect.

```julia
faces     = PolyVolume2D{Float64}[]
divisions = Tuple{Int,Int}[]

for j in 1:N_layers
    y_bot = layer_edges_norm[j]
    y_top = layer_edges_norm[j + 1]
    y_mid = (y_bot + y_top) / 2

    verts = SVector(
        Point2(0.0, y_bot), Point2(width, y_bot),
        Point2(width, y_top), Point2(0.0, y_top)
    )
    solidwalls = SVector((j == 1), true, false, true)

    # Exponential density decay
    ρ = exp(-y_mid / normalized_scale_height)

    # Spectral absorption: sigmoid from visible-transparent to IR-opaque
    sigmoid = [(1 / (1 + (4e-6 / λ_center[b])^6)) for b in 1:n_bins]
    layer_κ = [ρ * (κ_ir * sigmoid[b] + κ_vis * (1 - sigmoid[b]))
               for b in 1:n_bins]

    face = PolyVolume2D{Float64}(verts, solidwalls, n_bins, 1.0, 0.0)
    face.kappa_g   = layer_κ
    face.sigma_s_g = fill(0.0, n_bins)
    face.epsilon   = [fill(1.0, n_bins) for _ in 1:4]
    face.T_in_g    = -1.0       # solve for gas temperature
    face.q_in_g    = 0.0        # radiative equilibrium

    if j == 1
        face.T_in_w = [-1.0, 0.0, 0.0, 0.0]  # free surface, cold sides
    else
        face.T_in_w = [0.0, 0.0, 0.0, 0.0]   # cold boundaries
    end
    face.q_in_w = [0.0, 0.0, 0.0, 0.0]

    push!(faces, face)
    push!(divisions, (1, 2))
end
```

### Step 4: Add the solar source and build the mesh

```julia
sun_h_norm = sun_layer_height / L
verts_sun  = SVector(
    Point2(0.0, 1.0), Point2(width, 1.0),
    Point2(width, 1.0 + sun_h_norm), Point2(0.0, 1.0 + sun_h_norm)
)

face_sun = PolyVolume2D{Float64}(
    verts_sun, SVector(false, true, true, true), n_bins, κ_sun, 0.0);
face_sun.kappa_g   = fill(κ_sun, n_bins)
face_sun.sigma_s_g = fill(0.0, n_bins)
face_sun.epsilon   = [fill(1.0, n_bins) for _ in 1:4]
face_sun.T_in_g    = T_sun      # prescribed solar temperature
face_sun.q_in_g    = 0.0
face_sun.T_in_w    = [0.0, 0.0, 0.0, 0.0]
face_sun.q_in_w    = [0.0, 0.0, 0.0, 0.0]

push!(faces, face_sun);
push!(divisions, (1, 2));

mesh = RayTracingDomain2D(faces, divisions);
mesh.wavelength_band_limits = λ_edges
```

### Step 5: Ray trace and solve

```julia
mesh(2_000_000; method = :exchange)
solveEquilibrium!(mesh, mesh.F_smooth;
    max_iterations = 10_000, convergence_tol = 1e-12)
```

Ray tracing is performed independently for each spectral bin, computing separate exchange factor matrices that reflect the wavelength-dependent extinction. The spectral equilibrium solver then iterates to find the temperature distribution that simultaneously satisfies energy conservation across all bins.

### Step 6: Extract and plot the temperature profile

```julia
using Plots

gas_temps = Float64[]
altitudes = Float64[]

push!(gas_temps, mesh.fine_mesh[1][1].T_w[1])   # surface temperature
push!(altitudes, 0.0)

for j in 1:N_layers
    for k in 1:2
        push!(gas_temps, mesh.fine_mesh[j][k].T_g)
        push!(altitudes, mesh.fine_mesh[j][k].midPoint[2] * L)
    end
end

p = Plots.plot(gas_temps, altitudes ./ 1000,
    linewidth = 2, color = :black, marker = :circle, markersize = 3,
    ylabel = "Altitude (km)", xlabel = "Temperature (K)",
    title = "Atmospheric temperature profile\n(spectral greenhouse effect)",
    legend = false, dpi = 500,
    guidefontsize = 12, tickfontsize = 10,
    left_margin = 5Plots.mm, bottom_margin = 10Plots.mm,
    right_margin = 15Plots.mm)

savefig(p, "fig/spectral_greenhouse.png")
display(p)
```

![Spectral greenhouse atmosphere](fig/spectral_greenhouse.png)

The surface temperature emerges well above the bare blackbody equilibrium — a direct consequence of the spectral asymmetry between incoming (visible) and outgoing (infrared) radiation. Temperature decreases monotonically with altitude as the atmosphere thins and becomes transparent.

> **Note:** This is a simplified radiative equilibrium model without convection, latent heat, or detailed molecular absorption bands. Nevertheless, it captures the essential greenhouse mechanism from first principles: the spectral solver enforces energy conservation across the full spectrum to machine precision, and the temperature profile emerges purely from the exchange of radiation between layers.

The spectral solver is an unpublished extension of the grey GERT method described in [Bielefeld (2025)](https://arxiv.org/abs/2512.22157). It solves the coupled spectral equilibrium by iterating over Planck-weighted band contributions while preserving the exchange factor framework and its energy conservation guarantees.

---

## Example 3 — 3D Surface Enclosure

This example solves radiative equilibrium in a unit cube with transparent (non-participating) media. Two opposing faces have prescribed temperatures (1000 K and 0 K); the four side walls are in radiative equilibrium (unknown temperature, zero net heat flux). All surfaces are black (ε = 1). View factors are computed analytically using the method of Narayanaswamy (2015) — no ray tracing is needed.

### Step 1: Define the cube geometry

```julia
using RayTraceHeatTransfer
using GLMakie

# Cube vertices
points = [
    0.0 0.0 0.0;  # 1
    0.0 0.0 1.0;  # 2
    0.0 1.0 0.0;  # 3
    0.0 1.0 1.0;  # 4
    1.0 0.0 0.0;  # 5
    1.0 0.0 1.0;  # 6
    1.0 1.0 0.0;  # 7
    1.0 1.0 1.0   # 8
]

# Six faces (right-hand rule for outward normals)
faces = [
    1 2 4 3;  # Face 1 (x = 0) — hot
    5 6 8 7;  # Face 2 (x = 1) — cold
    1 5 7 3;  # Face 3 (z = 0)
    2 6 8 4;  # Face 4 (z = 1)
    3 4 8 7;  # Face 5 (y = 1)
    1 2 6 5   # Face 6 (y = 0)
]
```

### Step 2: Set boundary conditions and build the domain

```julia
Ndim = 11  # 11 × 11 subdivisions per face

epsilon = ones(size(faces, 1))                     # black surfaces
q_in_w  = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]           # zero net heat flux on sides
T_in_w  = [1000.0, 0.0, -1.0, -1.0, -1.0, -1.0]    # hot, cold, four unknown

domain3D = ViewFactorDomain3D(points, faces, Ndim, q_in_w, T_in_w, epsilon);
```

Faces 1 and 2 are the hot and cold walls at opposing ends of the cube. The four side faces have `T_in_w = -1.0` (unknown) and `q_in_w = 0.0` (radiative equilibrium), so their temperature distributions emerge from the solution.

### Step 3: Solve and visualise

```julia
solveEquilibrium!(domain3D, domain3D.F_smooth)

fig = Figure(size = (800, 700))
ax  = LScene(fig[1, 1], scenekw = (camera = cam3d!, show_axis = true))
plotField(ax, domain3D; field = :T)
fig
save("fig/3d_cube_temperature.png", fig, px_per_unit = 3)
```

![3D cube temperature field](fig/3d_cube_temperature.png)

The temperature field shows a smooth gradient from the hot face (1000 K) to the cold face (0 K), with the side walls at intermediate temperatures determined by radiative equilibrium. The analytical view factors ensure exact geometric accuracy without statistical noise.

### Reference

> Narayanaswamy, A. (2015). "An analytic expression for radiation view factor between two arbitrarily oriented planar polygons." *International Journal of Heat and Mass Transfer*, 91, 841–847.

---

## Documentation

The documentation of this package will gradually be rolled out in an online book format [here](https://gert.net/).

## References

The core methodology is presented in [Bielefeld (2025)](https://arxiv.org/abs/2512.22157). The 3D view factor implementation follows [Narayanaswamy (2015)](https://doi.org/10.1016/j.ijheatmasstransfer.2015.07.131). This work was inspired in part by [Howell, Mengüç, Daun & Siegel (2020)](https://www.routledge.com/Thermal-Radiation-Heat-Transfer/Howell-Menguc-Daun-Siegel/p/book/9780367347079).

## Authors

The primary author, developer and maintainer of this repository is Nikolaj Maack Bielefeld.

The functions for calculating 3D view factors analytically were originally written for MATLAB by Jacob A. Kerkhoff and Michael J. Wagner of University of Wisconsin-Madison, Energy Systems Optimization Lab, as described in [Kerkhoff & Wagner (2021)](https://asmedigitalcollection.asme.org/ES/proceedings-abstract/ES2021/84881/1114915).

## Declaration of AI Assistance

Parts of this package and its documentation were developed with assistance from Claude (Anthropic). All code, methods, and scientific content have been verified and validated by the author.