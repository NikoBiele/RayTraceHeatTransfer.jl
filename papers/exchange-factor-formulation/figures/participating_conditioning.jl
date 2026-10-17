# Figure: conditioning of the mixed-boundary system. A unit square with black walls at prescribed temperatures and
# every gas cell in radiative equilibrium, so all gas rows of M come from C = I − Fᵀ. With black walls and g = 1 on
# the gas rows, M does not depend on the single scattering albedo, so a non-scattering medium (κ = β) is used.
# The 1-norm condition number of M is plotted against the number of cells per side, for four optical thicknesses.

using RayTraceHeatTransfer                         # the GERT solver
using GeometryBasics, StaticArrays                 # points and fixed-size vectors for the geometry
using LinearAlgebra                                # cond, the 1-norm condition number
using CairoMakie                                   # Makie backend with vector (PDF) output
include("figure_style.jl")                         # shared scale, fonts and LaTeX strings for all figures of the paper
set_theme!(paper_makie_theme())                    # apply the paper's Makie theme to every figure drawn below

cells_per_side   = [5, 10, 20, 40]                 # mesh resolutions: N × N gas cells
optical_sides    = [0.1, 1.0, 10.0, 100.0]         # optical thickness of a side, τ_L = β · 1 m
rays_per_emitter = 2^12                            # rays per element, the same statistical quality on every mesh

vertices   = SVector(Point2(0.0, 0.0), Point2(1.0, 0.0), Point2(1.0, 1.0), Point2(0.0, 1.0))   # unit square [m]
solidWalls = SVector(true, true, true, true)       # all four walls are solid

conditioning = zeros(length(optical_sides), length(cells_per_side))   # cond₁(M) for every optical thickness (rows) and mesh (columns)

for (a, tau) in enumerate(optical_sides), (b, N) in enumerate(cells_per_side)
    face = PolyVolume2D{Float64}(vertices, solidWalls, 1, tau, 0.0)   # κ = β = τ_L on a 1 m side, σₛ = 0
    face.T_in_w  = [1000.0, 0.0, 0.0, 0.0]         # wall temperatures prescribed: wall rows of M come from D
    face.epsilon = [1.0, 1.0, 1.0, 1.0]            # black walls (b = 0)
    face.T_in_g  = -1.0                            # gas temperature unknown ...
    face.q_in_g  = 0.0                             # ... in radiative equilibrium: gas rows of M come from C

    mesh = RayTracingDomain2D([face], [(N, N)]; verbose = false)    # N × N cells
    emitters = length(mesh.surface_mapping) + length(mesh.volume_mapping)   # number of elements, 4N + N²
    println("τ_L = $tau, $N × $N cells ($emitters elements): tracing ...")  # progress, before the expensive step
    mesh(emitters * rays_per_emitter; method = :exchange, sampler = :sobol, verbose = false)   # exchange factors by ray tracing

    M = Matrix(buildSystemMatrix(mesh, mesh.F_raw))   # the mixed-boundary matrix exactly as the solver assembles it, dense
    conditioning[a, b] = cond(M, 1)                # ‖M‖₁ ‖M⁻¹‖₁
    println("τ_L = $tau, $N × $N cells: cond₁(M) = ", conditioning[a, b])   # reported in the text
end

markers    = ['+', '×', :utriangle, :rect]         # one marker per optical thickness: thin glyphs and hollow shapes
fills      = [:black, :black, :white, :white]      # glyphs are drawn in their colour; triangle and square are hollow
sizes      = [MARKER_SIZE, MARKER_SIZE, 14, 14]    # glyph size is a font size; the shapes need a smaller size to match
tau_labels = ["0.1", "1", "10", "100"]             # optical thicknesses as printed in the legend

fig = Figure(size = (SINGLE_WIDTH, 520))           # single-panel canvas, included at 0.6\linewidth
ax  = Axis(fig[1, 1];                              # the plot area
    xscale = log10, yscale = log10,                # both axes logarithmic
    xticks = (cells_per_side, string.(cells_per_side)),   # tick at every mesh
    yticks = (10.0 .^ (0:4), [L"10^{%$k}" for k in 0:4]), # tick at every decade
    xlabel = L"\mathrm{Cells\ per\ side}",
    ylabel = rich("cond", subscript("1"), "(", rich("M"; font = :bold), ")"))   # cond₁(M) with M upright and bold, as matrices are set in the paper
for (a, tau) in enumerate(optical_sides)           # one line per optical thickness
    scatterlines!(ax, cells_per_side, conditioning[a, :];
        color = :black, linewidth = 2,             # solid black line ...
        marker = markers[a], markersize = sizes[a],   # ... distinguished by marker
        markercolor = fills[a], strokecolor = :black, strokewidth = 1,   # thin outline on the hollow shapes
        label = L"\tau_L = %$(tau_labels[a])")     # legend entry
end
Legend(fig[0, 1], ax;                              # legend above the plot area, clear of the data in every figure
    orientation = :horizontal,                     # all entries in one row
    framevisible = false)                          # no frame: it sits outside the axes
fig                                                # show the figure for inspection

## Save, once the figure has been inspected
save("papers/exchange-factor-formulation/figures/participating_conditioning.pdf", fig)                    # vector output, named after the script
save("papers/exchange-factor-formulation/figures/participating_conditioning.png", fig; px_per_unit = 2)   # raster copy for inspection, at double resolution