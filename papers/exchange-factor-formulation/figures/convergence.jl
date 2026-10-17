# Figure: convergence in rays. The Example 1 problem (unit square, κ = 1 m⁻¹, bottom wall at 1000 K, all walls
# black, gas in radiative equilibrium) on an 11 × 11 mesh, whose centreline source function is compared to
# Crosbie & Schrenker (1984), as in Example 10 on gert.net but with the raw exchange factors (no smoothing).
# Every point is the mean over independent seeds, with error bars of one standard deviation across the seeds.

using RayTraceHeatTransfer                         # the GERT solver
using GeometryBasics, StaticArrays                 # points and fixed-size vectors for the geometry
using ConvolutionInterpolations                    # smooth interpolation of the tabulated reference
using Statistics                                   # mean and standard deviation across seeds
using CairoMakie                                   # Makie backend with vector (PDF) output
include("figure_style.jl")                         # shared scale, fonts and LaTeX strings for all figures of the paper
set_theme!(paper_makie_theme())                    # apply the paper's Makie theme to every figure drawn below

seeds = 1:10                                       # independent replicates: digital shifts for Sobol, streams for pseudorandom

function build_mesh(Ndim)                          # the Example 1 geometry and boundary conditions
    vertices = SVector(Point2(0.0, 0.0), Point2(1.0, 0.0), Point2(1.0, 1.0), Point2(0.0, 1.0))   # unit square
    face = PolyVolume2D{Float64}(vertices, SVector(true, true, true, true), 1, 1.0, 0.0)   # κ = 1, σₛ = 0
    face.T_in_w  = [1000.0, 0.0, 0.0, 0.0]        # bottom hot, rest cold
    face.epsilon = [1.0, 1.0, 1.0, 1.0]            # black walls
    face.T_in_g  = -1.0                            # unknown gas temperature
    face.q_in_g  = 0.0                             # radiative equilibrium
    return RayTracingDomain2D([face], [(Ndim, Ndim)]; verbose = false)
end

# Crosbie & Schrenker (1984): optical depth from the hot wall, and the source function S there
tau_ref = [0.0, 0.00611, 0.02037, 0.04251, 0.07216, 0.10884, 0.15194,
           0.20076, 0.25449, 0.31225, 0.37309, 0.43602, 0.50000, 0.56398,
           0.62691, 0.68775, 0.74551, 0.79924, 0.84806, 0.89116, 0.92784,
           0.95749, 0.97963, 0.99390, 1.00000]
S_ref   = [0.6293, 0.6198, 0.6017, 0.5767, 0.5460, 0.5108, 0.4724,
           0.4323, 0.3919, 0.3525, 0.3153, 0.2810, 0.2500, 0.2224,
           0.1981, 0.1768, 0.1584, 0.1424, 0.1287, 0.1171, 0.1073,
           0.0992, 0.0930, 0.0885, 0.0863]
S_ref_itp = convolution_interpolation(tau_ref, S_ref; kernel = :b13)   # the reference as a function of optical depth

# Average of the reference over the row between optical depths a and b, by the midpoint rule
function row_average(a, b; n = 1000)
    width = (b - a) / n                            # width of one subinterval
    total = sum(S_ref_itp(a + (k - 0.5) * width) for k in 1:n)   # the reference at every subinterval midpoint
    return total / n                               # its mean over the row
end

# Rms deviation of the centreline source function from Crosbie & Schrenker, for the Example 1
# problem on an Ndim × Ndim mesh, traced with `rays` rays by the given tracer, sampler and seed
function centreline_deviation(Ndim, rays; method = :exchange, sampler = :sobol, seed = 1)
    mesh = build_mesh(Ndim)                        # a fresh domain
    mesh(rays; method = method, sampler = sampler, seeds = seed, verbose = false)   # exchange factors by ray tracing, one replicate
    solveEquilibrium!(mesh, mesh.F_raw; verbose = false)   # single linear solve with the raw exchange factors (no smoothing in this paper)

    T = reshape([cell.T_g for cell in mesh.fine_mesh[1]], Ndim, Ndim)   # gas temperatures, first index along x
    S = (T ./ 1000.0) .^ 4                         # source function (T/T_hot)⁴
    c = div(Ndim + 1, 2)                           # the middle column

    # A cell holds the average over its width, which lies below the centreline value (the profile peaks there) by h²/24
    # times the magnitude of the curvature across the centreline; the second difference of the three middle columns
    # estimates it, and the offset is added back
    S_centre = S[c, :] .- (S[c-1, :] .- 2 .* S[c, :] .+ S[c+1, :]) ./ 24

    edges = range(0, 1, length = Ndim + 1)         # row boundaries in optical depth
    S_ref_rows = [row_average(edges[j], edges[j+1]) for j in 1:Ndim]   # the reference averaged over each row
    return sqrt(sum(abs2, S_centre .- S_ref_rows) / Ndim)              # rms deviation
end

# The ways of computing the exchange factors: tracer, sampler, and a label for the figure
configurations = [(:exchange, :random, L"\mathrm{Pseudorandom}"),
                  (:exchange, :sobol,  L"\mathrm{Sobol}")]

Ndim_rays        = 11                              # mesh of the ray study
emitters         = 4Ndim_rays + Ndim_rays^2        # number of elements on that mesh
rays_per_emitter = [2^10, 2^12, 2^14, 2^16, 2^18, 2^20] # rays per emitter

# rms deviation for every configuration, number of rays per emitter and seed
ray_study = [centreline_deviation(Ndim_rays, emitters * r; method = m, sampler = s, seed = k)
             for (m, s, _) in configurations, r in rays_per_emitter, k in seeds]

ray_mean = dropdims(mean(ray_study, dims = 3), dims = 3)   # mean over seeds, per configuration and ray count
ray_std  = dropdims(std(ray_study, dims = 3), dims = 3)    # standard deviation over seeds

println("rays, mean over seeds: ", ray_mean)       # reported in the text
println("rays, standard deviation over seeds: ", ray_std)

markers = ['×', '+']                               # thin × for pseudorandom, thin + for Sobol

fig = Figure(size = (SINGLE_WIDTH, 520))           # single-panel canvas, included at 0.6\linewidth
ax  = Axis(fig[1, 1];                              # the plot area
    xscale = log10, yscale = log10,                # both axes logarithmic
    xticks = (rays_per_emitter, [L"2^{%$k}" for k in 10:2:20]),   # tick at every ray count
    xlabel = L"\mathrm{Rays\ per\ emitter}",
    ylabel = L"\mathrm{rms\ deviation\ of\ } S")
for (k, (_, _, label)) in enumerate(configurations)
    low = min.(ray_std[k, :], 0.9 .* ray_mean[k, :])   # lower error bar, capped so it stays positive on the log axis
    errorbars!(ax, rays_per_emitter, ray_mean[k, :], low, ray_std[k, :];
               color = :black, whiskerwidth = 10)      # one standard deviation across seeds
    scatterlines!(ax, rays_per_emitter, ray_mean[k, :]; color = :black, linewidth = 2,
                  marker = markers[k], markersize = MARKER_SIZE, label = label)   # mean over seeds, with its legend entry
end
guide = ray_mean[1, 1] .* (rays_per_emitter ./ rays_per_emitter[1]) .^ (-1 / 2)   # classical Monte Carlo rate, as in Example 10
lines!(ax, rays_per_emitter, guide; color = :black, linestyle = :dash, linewidth = 2,
       label = L"\propto \, N^{-1/2}")                # dashed guide line, with its legend entry
ylims!(ax, 1e-4, 1e-1)                             # y from 10⁻⁴ to 10⁻¹
ax.yticks = (10.0 .^ (-4:-1), [L"10^{%$k}" for k in -4:-1])   # a tick at every decade
Legend(fig[0, 1], ax;                              # legend above the plot area, clear of the data in every figure
    orientation = :horizontal,                     # all entries in one row
    framevisible = false)                          # no frame: it sits outside the axes
fig                                                # show the figure for inspection

## Save, once the figure has been inspected
save("papers/exchange-factor-formulation/figures/convergence.pdf", fig)                    # vector output, named after the script
save("papers/exchange-factor-formulation/figures/convergence.png", fig; px_per_unit = 2)   # raster copy for inspection, at double resolution