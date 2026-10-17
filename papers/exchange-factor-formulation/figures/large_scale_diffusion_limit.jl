# Figure: large-scale diffusion limit. A 1000 m × 1 m rectangle (approximately one-dimensional) filled with an
# absorbing-emitting, non-scattering medium (β = 100 m⁻¹), one long wall hot, all walls black, meshed with
# 1001 × 1001 cells. The exchange factors are traced with pseudorandom and with Sobol-sampled rays and used
# as traced; the centreline emissive power is compared to the diffusion approximation between parallel plates.
# The smoothed solutions are computed for reference only and are printed, not plotted.

using RayTraceHeatTransfer                         # the GERT solver
using GeometryBasics, StaticArrays                 # points and fixed-size vectors for the geometry
using SparseArrays                                 # nnz of the sparse exchange factor matrix
using CairoMakie                                   # Makie backend with vector (PDF) output
include("figure_style.jl")                         # shared scale, fonts and LaTeX strings for all figures of the paper
set_theme!(paper_makie_theme())                    # apply the paper's Makie theme to every figure drawn below

T_hot = 1000.0                                     # temperature of the hot wall [K]
beta  = 100.0                                      # extinction coefficient [m⁻¹], equal to κ since σₛ = 0
L     = 1.0                                        # distance between the two long walls [m]
Nx    = 1001                                       # cells along the 1000 m length
Ny    = 1001                                       # cells across the 1 m gap
rays  = (Nx*Ny+2*Nx+2*Ny)*2^10                     # total number of rays traced per sampler

# Normalized emissive power (T/T_hot)⁴ along the middle column of the solved mesh, from the hot wall upwards
function centreline_S(mesh, nx, ny)
    T_cells = reshape([cell.T_g for cell in mesh.fine_mesh[1]], nx, ny)    # gas temperatures, first index along x, second across the gap
    return (T_cells[div(nx + 1, 2), :] ./ T_hot) .^ 4                      # the middle column, normalized
end

# One complete run: build the mesh, trace the exchange factors with the given sampler, solve once with the raw
# exchange factors, then smooth and solve again, and return both centreline profiles with size, sparsity and timings
function run_case(sampler, nx, ny, nrays)
    vertices   = SVector(Point2(0.0, 0.0), Point2(1000.0, 0.0), Point2(1000.0, L), Point2(0.0, L))   # the rectangle, walls ordered bottom, right, top, left
    solidWalls = SVector(true, true, true, true)   # all four walls are solid
    face = PolyVolume2D{Float64}(vertices, solidWalls, 1, beta, 0.0)   # one spectral bin, κ = β, σₛ = 0
    face.T_in_w  = [T_hot, 0.0, 0.0, 0.0]          # bottom (long) wall hot, the others at 0 K
    face.epsilon = [1.0, 1.0, 1.0, 1.0]            # black walls
    face.T_in_g  = -1.0                            # gas temperature unknown ...
    face.q_in_g  = 0.0                             # ... in radiative equilibrium

    mesh = RayTracingDomain2D([face], [(nx, ny)])  # nx × ny cells
    t_trace = @elapsed mesh(nrays; method = :exchange, sampler = sampler)   # exchange factors by ray tracing, timed
    t_solve = @elapsed solveEquilibrium!(mesh, mesh.F_raw)                  # single linear solve with the raw exchange factors, timed
    S_raw   = centreline_S(mesh, nx, ny)           # centreline profile of the raw solution
    err_raw = mesh.energy_error                    # relative energy conservation error of the raw solution

    F = mesh.F_raw                                 # the exchange factor matrix as traced
    N   = size(F, 1)                               # number of elements, walls and gas
    nzF = F isa AbstractSparseMatrix ? nnz(F) : count(!iszero, F)          # stored non-zero exchange factors

    t_smooth = @elapsed smooth!(mesh)              # reciprocity and closure enforced by the package smoother, timed
    solveEquilibrium!(mesh, mesh.F_smooth)         # the same solve with the smoothed exchange factors
    S_smooth = centreline_S(mesh, nx, ny)          # centreline profile of the smoothed solution

    return (S_raw = S_raw, S_smooth = S_smooth,    # the two centreline profiles
            N = N, nnz = nzF,                      # size and sparsity of F
            t_trace = t_trace, t_solve = t_solve, t_smooth = t_smooth,   # timings [s]
            sol_j = mesh.J,
            energy_error = err_raw)                # relative energy conservation error of the raw solution
end

# Warm-up on a small mesh, so that the timed runs below exclude Julia's just-in-time compilation
run_case(:random, 3, 51, 10^6)                     # compiles the pseudorandom tracer, the solver and the smoother
run_case(:sobol,  3, 51, 10^6)                     # compiles the Sobol tracer

res_random = run_case(:random, Nx, Ny, rays)       # pseudorandom rays, i.e. classical Monte Carlo
GC.gc()                                            # release the first mesh before the second trace
res_sobol  = run_case(:sobol,  Nx, Ny, rays)       # Sobol-sampled rays
GC.gc()                                            # release the second mesh

# Diffusion approximation between infinite parallel plates with temperature-jump boundary conditions:
# the net flux q follows from the wall emissive powers, and the gas emissive power falls linearly across the gap
function diffusion_emissive_power(z, β, D, ε1, ε2, E1, E2)
    q  = (E1 - E2) / (3β * D / 4 + 1 / ε1 + 1 / ε2 - 1)   # net radiative flux between the plates [W/m²]
    E0 = E1 + q * (1 / 2 - 1 / ε1)                          # gas emissive power just off the hot wall (the jump)
    return E0 - (3β * z / 4) * q                            # linear profile across the gap
end

σ_SB    = 5.670374419e-8                           # Stefan–Boltzmann constant [W/(m² K⁴)]
E_hot   = σ_SB * T_hot^4                           # blackbody emissive power of the hot wall [W/m²]
z_cells = ((1:Ny) .- 0.5) .* (L / Ny)              # distance of each cell centre from the hot wall [m]
z_line  = range(0.0, L, length = 201)              # fine grid for the reference curve [m]
S_line  = diffusion_emissive_power.(z_line, beta, L, 1.0, 1.0, E_hot, 0.0) ./ E_hot    # reference curve, normalized
S_ref   = diffusion_emissive_power.(z_cells, beta, L, 1.0, 1.0, E_hot, 0.0) ./ E_hot   # reference at the cell centres
rms(S)  = sqrt(sum(abs2, S .- S_ref) / Ny)         # rms deviation of a centreline profile from the diffusion approximation

# Numbers reported in the text, for each sampler; the smoothed rms is for reference only
for (name, res) in (("pseudorandom", res_random), ("Sobol", res_sobol))
    println(name, ":")                                                # header for this sampler
    println("  elements N:                         ", res.N)          # matrix dimension
    println("  non-zero exchange factors:          ", res.nnz)        # stored entries of F
    println("  fraction of N² stored:              ", round(res.nnz / res.N^2; sigdigits = 2))   # sparsity of F
    println("  storage of a dense N × N matrix:    ", round(8 * res.N^2 / 1e12; sigdigits = 2), " TB")   # Float64, for comparison
    println("  storage of the dense gas-gas block: ", round(8 * (Nx * Ny)^2 / 1e12; sigdigits = 2), " TB")   # Noble's total exchange areas between gas elements, dense once the medium scatters
    println("  tracing time:                       ", round(res.t_trace; digits = 1), " s")   # time for the rays
    println("  solving time:                       ", round(res.t_solve; digits = 1), " s")   # single linear solve and post-processing
    println("  smoothing time (reference only):    ", round(res.t_smooth; digits = 1), " s")  # package smoother
    println("  rms deviation, raw F:               ", round(rms(res.S_raw); sigdigits = 2))      # reported in the text
    println("  rms deviation, smoothed F:          ", round(rms(res.S_smooth); sigdigits = 2))   # reference only
    println("  relative energy conservation error: ", res.energy_error)                      # machine precision expected
    println("  minimum j:                          ", minimum(res.sol_j))
end

fig = Figure(size = (TEXT_WIDTH, 520))     # two-panel canvas, included at \linewidth with the same font scale as the single panels
ax_random = Axis(fig[1, 1];                        # left panel: pseudorandom rays
    title  = L"\mathrm{Pseudorandom}",             # panel title
    xlabel = L"\tau/\tau_L",                       # optical depth from the hot wall, relative to the total
    ylabel = L"(T/T_\mathrm{hot})^4")              # normalized emissive power
ax_sobol = Axis(fig[1, 2];                         # right panel: Sobol-sampled rays
    title  = L"\mathrm{Sobol}",                    # panel title
    xlabel = L"\tau/\tau_L")                       # optical depth from the hot wall, relative to the total
linkyaxes!(ax_random, ax_sobol)                    # the same vertical scale in both panels
hideydecorations!(ax_sobol; grid = false, ticks = false)   # the shared axis is labelled once, on the left

for (ax, res) in ((ax_random, res_random), (ax_sobol, res_sobol))
    scatter!(ax, z_cells ./ L, res.S_raw;          # the solution with raw exchange factors at every cell centre, drawn first
        marker = :circle,                          # circles ...
        color = :transparent,                      # ... hollow, so that overlapping markers stay distinguishable
        strokecolor = :black,                      # thin black outline
        strokewidth = 0.5,                         # outline width
        markersize = MARKER_SIZE,                  # marker size, the same in every Makie figure
        label = L"\mathrm{Proposed\ formulation}") # legend entry
    lines!(ax, z_line ./ L, S_line;                # the diffusion approximation, first as a wider black line ...
        color = :black, linewidth = 4,             # ... which shows as a thin outline around the white line drawn next
        label = L"\mathrm{Diffusion\ approximation}")   # same label as the white line, merged into one legend entry
    lines!(ax, z_line ./ L, S_line;                # ... then as a white line on top, cutting through the scatter
        color = :white, linewidth = 2,             # white core, visible against the dense black markers
        label = L"\mathrm{Diffusion\ approximation}")   # same label as the outline, merged into one legend entry
end
Legend(fig[0, 1:2], ax_random;                     # one legend above both panels, clear of the data
    orientation = :horizontal,                     # all entries in one row
    framevisible = false,                          # no frame: it sits outside the axes
    merge = true)                                  # entries with the same label are overlaid, so the line shows outlined
fig                                                # show the figure for inspection

## Save, once the figure has been inspected
save("papers/exchange-factor-formulation/figures/large_scale_diffusion_limit.pdf", fig)                    # vector output, named after the script
save("papers/exchange-factor-formulation/figures/large_scale_diffusion_limit.png", fig; px_per_unit = 3)   # raster copy for inspection, at double resolution

############################################
############################################

## Graphical abstract in Elsevier's 5:2 frame: the solution on the left, the storage it needed on the right
GA_FONT   = 36                                                       # one font size for every text element, large enough for the 200-pixel-high display
dense_GB  = 8 * (Nx * Ny)^2 / 1e9                                    # gas-gas block of Noble's total exchange areas, dense Float64 [GB], strictly positive once the medium scatters
sparse_GB = (16 * res_sobol.nnz + 8 * (res_sobol.N + 1)) / 1e9      # SparseMatrixCSC{Float64,Int64}: values + row indices + column pointers [GB]
ratio     = Int(round(dense_GB / sparse_GB; sigdigits = 2))          # how many times smaller the sparse storage is

fig_ga = Figure(size = (TEXT_WIDTH, round(Int, TEXT_WIDTH / 2.5)))   # Elsevier's 5:2 graphical abstract frame

ax_sol = Axis(fig_ga[1, 1];                                          # left: the solution against the reference
    title  = L"N \approx 10^6\ \mathrm{elements}",                   # problem size
    xlabel = L"\tau/\tau_L",                                         # optical depth from the hot wall, relative to the total
    ylabel = L"(T/T_\mathrm{hot})^4",                                # normalized emissive power
    titlesize = GA_FONT,                                             # title size
    xlabelsize = GA_FONT, ylabelsize = GA_FONT,                      # axis label sizes
    xticklabelsize = GA_FONT, yticklabelsize = GA_FONT)              # tick label sizes
scatter!(ax_sol, z_cells ./ L, res_sobol.S_raw;                      # the solution with raw exchange factors at every cell centre, drawn first
    marker = :circle,                                                # circles ...
    color = :transparent,                                            # ... hollow, so that overlapping markers stay distinguishable
    strokecolor = :black,                                            # thin black outline
    strokewidth = 0.5,                                               # outline width
    markersize = MARKER_SIZE,                                        # marker size, the same in every Makie figure
    label = L"\mathrm{Proposed\ formulation}")                       # legend entry
lines!(ax_sol, z_line ./ L, S_line;                                  # the diffusion approximation, first as a wider black line ...
    color = :black, linewidth = 4,                                   # ... which shows as a thin outline around the white line drawn next
    label = L"\mathrm{Diffusion\ approximation}")                    # same label as the white line, merged into one legend entry
lines!(ax_sol, z_line ./ L, S_line;                                  # ... then as a white line on top, cutting through the scatter
    color = :white, linewidth = 2,                                   # white core, visible against the dense black markers
    label = L"\mathrm{Diffusion\ approximation}")                    # same label as the outline, merged into one legend entry
axislegend(ax_sol;                                                   # legend in the empty upper-right corner of the plot
    position = :rt,                                                  # right, top
    framevisible = false,                                            # no frame
    labelsize = GA_FONT-4,                                             # legend text size
    merge = true)                                                    # the two line entries are overlaid into one outlined line

ax_mem = Axis(fig_ga[1, 2];                                          # right: storage of the dense total exchange areas against the sparse exchange factor matrix
    title  = "Storage: $(ratio)× less",                               # the message of the panel: the ratio of the two bars
    titlefont = :regular,                                            # axis titles default to bold; keep the title upright like the rest of the panel
    ylabel = L"\mathrm{TB}",                                         # unit, as in the text of the paper
    xticks = ([1, 2], [L"\overline{\textbf{GG}}", L"\textbf{F}"]),   # the dense gas-gas block of the total exchange areas, in the notation of the paper, against the sparse exchange factor matrix; \textbf for upright bold, since \mathbf is sloped in MathTeXEngine
    limits = (0.4, 2.6, 0, 1.25 * dense_GB / 1000),                  # linear axis from zero, with room above the dense bar for its label
    titlesize = GA_FONT,                                             # title size
    ylabelsize = GA_FONT,                                            # axis label size
    xticklabelsize = GA_FONT, yticklabelsize = GA_FONT)              # tick label sizes
barplot!(ax_mem, [1, 2], [dense_GB, sparse_GB] ./ 1000;              # one bar per storage format, in TB
    color = [:gray70, :black],                                       # dense in grey, sparse in black
    width = 0.6)                                                     # bar width
text!(ax_mem, 1, dense_GB / 1000;                                    # label above the dense bar
    text = L"%$(round(dense_GB / 1000; digits = 1))\ \mathrm{TB}",   # in TB
    fontsize = GA_FONT,                                              # label size
    align = (:center, :bottom), offset = (0, 6))                     # centred just above the bar
text!(ax_mem, 2, sparse_GB / 1000;                                   # label where the sparse bar is, just above the axis
    text = L"%$(round(sparse_GB; digits = 1))\ \mathrm{GB}",         # in GB, since in TB it would round to zero
    fontsize = GA_FONT,                                              # label size
    align = (:center, :bottom), offset = (0, 6))                     # centred just above the (invisible) bar
colsize!(fig_ga.layout, 1, Relative(0.6))                            # the solution gets the larger share of the width
fig_ga                                                               # show the figure for inspection

## Save the graphical abstract, once it has been inspected
save("papers/exchange-factor-formulation/figures/large_scale_graphical_abstract.pdf", fig_ga)                    # vector output
save("papers/exchange-factor-formulation/figures/large_scale_graphical_abstract.png", fig_ga; px_per_unit = 3)   # raster copy, 3600 × 1440 pixels, above Elsevier's minimum