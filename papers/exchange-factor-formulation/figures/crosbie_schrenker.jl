# Figure: comparison with Crosbie & Schrenker (1984). Square media of side 1 m, exposed to diffuse radiation from a
# black wall at T_hot, the other walls black and cold, the medium scattering isotropically and its emission neglected
# (gas prescribed at 0 K). (a) albedo series at τ_L = 1; (b) optical-thickness series at ω = 1. The quantity compared
# is the centreline mean incident intensity, S/ω = g/(w σT_hot⁴), the incident power per unit emission capacity.

using RayTraceHeatTransfer                         # the GERT solver
using GeometryBasics, StaticArrays                 # points and fixed-size vectors for the geometry
using ConvolutionInterpolations                    # smooth interpolation of the tabulated reference
using CairoMakie                                   # Makie backend with vector (PDF) output
include("figure_style.jl")                         # shared scale, fonts and LaTeX strings for all figures of the paper
set_theme!(paper_makie_theme())                    # apply the paper's Makie theme to every figure drawn below

T_hot = 1000.0                                     # temperature of the illuminating wall [K]
N     = 21                                         # cells per side
σ_SB  = 5.670374419e-8                             # Stefan–Boltzmann constant [W/(m² K⁴)]
rays_per_emitter(tau) = round(Int, 2^16 * max(1.0, 1.0 / tau))   # rays per element; thin media need more, since a ray
                                                                  # interacts in a thin cell with probability ∝ its optical thickness

# Reference: Crosbie, A. L. and Schrenker, R. G. (1984), Radiative transfer in a two-dimensional rectangular medium
# exposed to diffuse radiation, JQSRT, doi:10.1016/0022-4073(84)90095-5. Centreline source function S/ω against the
# optical depth from the illuminated wall, relative to the total (τ/τ_L), for the six cases compared here.
tau_ratio_thin  = [0.0, 0.00611, 0.02037, 0.04251, 0.07216, 0.10884, 0.15194, 0.20076, 0.25449, 0.31225, 0.37309,
                   0.43602, 0.50000, 0.56398, 0.62691, 0.68775, 0.74551, 0.79924, 0.84806, 0.89116, 0.92784, 0.95749,
                   0.97963, 0.99390, 1.00000]      # optical depths of the tables for τ_L ≤ 1
tau_ratio_thick = [0.0, 0.00611, 0.02037, 0.04251, 0.07216, 0.10884, 0.15194, 0.20076, 0.25449, 0.31225, 0.37309,
                   0.43602, 0.50000, 0.56398, 0.62691, 0.68775, 0.74551, 0.79924, 0.84806, 0.89116, 0.92784, 0.95749,
                   0.97963, 0.99389, 1.00000]      # optical depths of the tables for τ_L = 4 and 10
reference = Dict(                                  # (τ_L, ω) => (τ/τ_L, S/ω)
    (0.25, 1.0) => (tau_ratio_thin,  [0.5415, 0.5358, 0.5238, 0.5062, 0.4837, 0.4571, 0.4274, 0.3958, 0.3636, 0.3321,
                                      0.3023, 0.2748, 0.2500, 0.2280, 0.2088, 0.1921, 0.1777, 0.1654, 0.1550, 0.1463,
                                      0.1392, 0.1335, 0.1292, 0.1263, 0.1249]),
    (1.0, 1.0)  => (tau_ratio_thin,  [0.6293, 0.6198, 0.6017, 0.5767, 0.5460, 0.5108, 0.4724, 0.4323, 0.3919, 0.3525,
                                      0.3153, 0.2810, 0.2500, 0.2224, 0.1981, 0.1768, 0.1584, 0.1424, 0.1287, 0.1171,
                                      0.1073, 0.0992, 0.0930, 0.0885, 0.0863]),
    (1.0, 0.9)  => (tau_ratio_thin,  [0.6102, 0.5996, 0.5798, 0.5530, 0.5208, 0.4844, 0.4452, 0.4048, 0.3646, 0.3258,
                                      0.2896, 0.2566, 0.2270, 0.2009, 0.1782, 0.1585, 0.1416, 0.1271, 0.1148, 0.1044,
                                      0.0958, 0.0888, 0.0833, 0.0795, 0.0776]),
    (1.0, 0.5)  => (tau_ratio_thin,  [0.5510, 0.5367, 0.5117, 0.4797, 0.4431, 0.4037, 0.3632, 0.3231, 0.2846, 0.2488,
                                      0.2164, 0.1878, 0.1631, 0.1419, 0.1240, 0.1090, 0.0965, 0.0862, 0.0777, 0.0707,
                                      0.0652, 0.0608, 0.0575, 0.0553, 0.0543]),
    (4.0, 1.0)  => (tau_ratio_thick, [0.7948, 0.7799, 0.7526, 0.7155, 0.6706, 0.6196, 0.5642, 0.5066, 0.4489, 0.3931,
                                      0.3407, 0.2928, 0.2500, 0.2124, 0.1798, 0.1517, 0.1278, 0.1074, 0.0902, 0.0757,
                                      0.0636, 0.0538, 0.0461, 0.0407, 0.0378]),
    (10.0, 1.0) => (tau_ratio_thick, [0.8954, 0.8801, 0.8511, 0.8099, 0.7582, 0.6979, 0.6313, 0.5610, 0.4903, 0.4219,
                                      0.3582, 0.3006, 0.2500, 0.2063, 0.1691, 0.1377, 0.1114, 0.0894, 0.0711, 0.0559,
                                      0.0435, 0.0335, 0.0259, 0.0207, 0.0180]),
)
optical_sides = [0.25, 1.0, 4.0, 10.0]             # τ_L of the traced geometries, β = τ_L on a 1 m side
albedos_at    = Dict(0.25 => [1.0], 1.0 => [1.0, 0.9, 0.5], 4.0 => [1.0], 10.0 => [1.0])   # cases solved on each geometry

vertices   = SVector(Point2(0.0, 0.0), Point2(1.0, 0.0), Point2(1.0, 1.0), Point2(0.0, 1.0))   # unit square, walls ordered bottom, right, top, left
solidWalls = SVector(true, true, true, true)       # all four walls are solid

# A cell holds the average over its width, which lies below the centreline value (the profile peaks there) by h²/24
# times the magnitude of the curvature across the centreline; the second difference of the three middle columns
# estimates it, and the offset is added back
function centreline_incident(m, beta)
    S_over_omega = [cell.g_g / (4 * beta * cell.volume * σ_SB * T_hot^4) for cell in m.fine_mesh[1]]   # every cell, row by row from the bottom
    S_over_omega = reshape(S_over_omega, N, N)                           # first index: column (x), second index: row (y)
    c = div(N + 1, 2)                              # the middle column
    return S_over_omega[c, :] .- (S_over_omega[c-1, :] .- 2 .* S_over_omega[c, :] .+ S_over_omega[c+1, :]) ./ 24   # centreline values, from the illuminated wall upwards
end

# The reference averaged over each row of cells, by the midpoint rule on an interpolation of the table
function row_averaged(tau_ref, S_over_omega_ref; n = 1000)
    itp   = convolution_interpolation(tau_ref, S_over_omega_ref; kernel = :b13)   # the reference as a function of τ/τ_L
    edges = range(0.0, 1.0, length = N + 1)        # row boundaries in τ/τ_L
    return [sum(itp(edges[k] + (i - 0.5) * (edges[k+1] - edges[k]) / n) for i in 1:n) / n for k in 1:N]   # mean over each row
end

computed = Dict{Tuple{Float64,Float64},Vector{Float64}}()   # (τ_L, ω) => computed centreline S/ω

for tau in optical_sides
    face = PolyVolume2D{Float64}(vertices, solidWalls, 1, tau, 0.0)   # β = τ_L on a 1 m side; the split is set per case
    face.T_in_w  = [T_hot, 0.0, 0.0, 0.0]          # bottom wall illuminates diffusely, the others are cold
    face.epsilon = [1.0, 1.0, 1.0, 1.0]            # black walls
    face.T_in_g  = 0.0                             # gas prescribed at 0 K: its emission is neglected, as in the reference

    mesh = RayTracingDomain2D([face], [(N, N)]; verbose = false)    # N × N cells
    emitters = length(mesh.surface_mapping) + length(mesh.volume_mapping)   # number of elements
    println("τ_L = $tau: tracing ...")             # progress, before the expensive step
    mesh(emitters * rays_per_emitter(tau); method = :exchange, sampler = :sobol, verbose = false)   # exchange factors, reused for every albedo

    for omega in albedos_at[tau]                   # every albedo compared on this geometry
        for cell in mesh.fine_mesh[1]              # same extinction, so the same exchange factors ...
            cell.kappa_g   = (1 - omega) * tau     # ... with absorption κ = (1 − ω)β ...
            cell.sigma_s_g = omega * tau           # ... and scattering σₛ = ωβ
        end
        solveEquilibrium!(mesh, mesh.F_raw; verbose = false)   # single linear solve with the raw exchange factors
        S_over_omega = centreline_incident(mesh, tau)         # computed centreline
        computed[(tau, omega)] = S_over_omega                 # keep it for the figure
        tau_ref, S_over_omega_ref = reference[(tau, omega)]   # the matching table
        rms = sqrt(sum(abs2, S_over_omega .- row_averaged(tau_ref, S_over_omega_ref)) / N)   # rms deviation from the row-averaged reference
        println("τ_L = $tau, ω = $omega: rms deviation = ", round(rms; sigdigits = 2),
                ", relative energy conservation error = ", mesh.energy_error)   # reported in the text
    end
end

tau_cells = ((1:N) .- 0.5) ./ N                    # τ/τ_L at the cell centres

markers      = ['+', '×', :utriangle, :rect]       # thin + and × glyphs, then hollow triangle and square
marker_fills = [:black, :black, :white, :white]    # glyphs in black; triangle and square filled white, i.e. hollow
marker_sizes = [MARKER_SIZE, MARKER_SIZE, 14, 14]  # glyphs at the shared size; the shapes smaller, to match them visually
stroke_widths = [0.0, 0.0, 1.5, 1.5]               # outline only on the hollow shapes (an outline would make a glyph bold)

fig = Figure(size = (TEXT_WIDTH, 600))             # full text width; height for two panels with their legends above

# One panel: the reference tables as solid lines, the computed centrelines as black markers, and a legend above
function panel!(column, cases, labels, title)
    ax = Axis(fig[1, column]; xlabel = L"\tau/\tau_L", ylabel = L"S/\omega", title = title,
              limits = (0, 1, 0, 1))               # common axes for both panels
    for (k, case) in enumerate(cases)
        tau_ref, S_over_omega_ref = reference[case]           # the table for this case
        lines!(ax, tau_ref, S_over_omega_ref; color = :black, linewidth = 2,
               label = k == 1 ? L"\mathrm{Crosbie\ and\ Schrenker}" : nothing)   # one legend entry for all reference lines
        scatter!(ax, tau_cells, computed[case];
                 marker = markers[k], color = marker_fills[k],          # marker shape and fill for this case
                 strokecolor = :black, strokewidth = stroke_widths[k],  # black outline on the hollow shapes only
                 markersize = marker_sizes[k], label = labels[k])       # legend entry for this case
    end
    Legend(fig[0, column], ax;                     # legend above the panel, clear of the data as in every figure
        orientation = :horizontal,                 # entries side by side ...
        nbanks = 2,                                # ... in two rows, since the panel is only half the text width
        framevisible = false)                      # no frame: it sits outside the axes
    return ax
end

panel!(1, [(1.0, 1.0), (1.0, 0.9), (1.0, 0.5)],    # (a) albedo series at τ_L = 1
       [L"\omega = 1", L"\omega = 0.9", L"\omega = 0.5"], L"\tau_L = 1")
panel!(2, [(0.25, 1.0), (1.0, 1.0), (4.0, 1.0), (10.0, 1.0)],   # (b) optical-thickness series at ω = 1
       [L"\tau_L = 0.25", L"\tau_L = 1", L"\tau_L = 4", L"\tau_L = 10"], L"\omega = 1")
fig                                                # show the figure for inspection

## Save, once the figure has been inspected
save("papers/exchange-factor-formulation/figures/crosbie_schrenker.pdf", fig)                    # vector output, named after the script
save("papers/exchange-factor-formulation/figures/crosbie_schrenker.png", fig; px_per_unit = 2)   # raster copy for inspection, at double resolution