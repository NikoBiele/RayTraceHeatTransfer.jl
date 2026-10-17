# Figure: symmetry validation on a triangulated sphere. Black icosphere enclosures at two refinement levels, with a
# hot cap (1000 K) at the top, an antipodal cold cap (0 K) at the bottom, and all other triangles in radiative
# equilibrium. Because the system is linear in emissive power and the sphere is centrally symmetric, swapping the
# caps maps every triangle onto its antipode, so T_k⁴ + T_{-k}⁴ = T_hot⁴ + T_cold⁴ exactly for every equilibrium
# triangle. Left: the temperature field on the finest sphere. Right: the largest deviation from the identity, in
# Float64 and in 128-bit arithmetic: both sit at their own machine precision, so the deviation is round-off.

using RayTraceHeatTransfer                         # the GERT solver
using CairoMakie                                   # Makie backend with vector (PDF) output
include("figure_style.jl")                         # shared scale, fonts and LaTeX strings for all figures of the paper
set_theme!(paper_makie_theme())                    # apply the paper's Makie theme to every figure drawn below
include(joinpath(pkgdir(RayTraceHeatTransfer), "examples", "icosphere_mesh.jl"))   # icosphere generator from the package examples

levels     = 1:2                                   # icosphere refinement levels: 80, 320 and 1280 triangles (level 0 has no antipodal caps)
levels_big = 1:2                                   # levels run in 128-bit arithmetic (1:2 for a quick look; level 3 is the costly one)
n_cap      = 6                                     # triangles per cap
BITS       = 128                                   # BigFloat precision of the second series [bits]

# Solve the cap problem on each level in the floating-point type T and return the number of triangles,
# the largest relative deviation from the symmetry identity, and the solved domain of the finest level
function symmetry_deviation(::Type{T}, levels) where {T}
    T_hot  = T(1000)                               # temperature of the hot cap [K]
    T_cold = T(0)                                  # temperature of the cold cap [K]
    n_triangles = Int[]                            # number of triangles at each level
    deviation   = Float64[]                        # largest relative violation of the identity at each level
    finest      = nothing                          # the solved domain of the finest level

    for level in levels
        points, faces = icosphere_mesh(T, level)   # the triangulated sphere, built in T
        n_tri = size(faces, 1)                     # number of triangles
        centroids = [vec(sum(points[faces[i, :], :], dims = 1)) ./ 3 for i in 1:n_tri]   # triangle centres
        z = [c[3] for c in centroids]              # their heights
        hot  = partialsortperm(z, 1:n_cap, rev = true)   # the highest triangles form the hot cap
        cold = partialsortperm(z, 1:n_cap)         # the lowest form the cold cap

        antipode = [argmin([sum(abs2, centroids[j] .+ centroids[i]) for j in 1:n_tri]) for i in 1:n_tri]   # each triangle's opposite triangle
        Set(antipode[hot]) == Set(cold) ||         # the identity holds only if the caps are exact antipodes
            error("level $level: the hot and cold caps are not antipodal; choose another cap size")

        T_in_w = fill(-one(T), n_tri)              # all temperatures unknown ...
        T_in_w[hot]  .= T_hot                      # ... except the hot cap ...
        T_in_w[cold] .= T_cold                     # ... and the cold cap
        q_in_w  = zeros(T, n_tri)                  # zero net heat flux on the unknown triangles: radiative equilibrium
        epsilon = ones(T, n_tri)                   # black surfaces

        println("$T, level $level: $n_tri triangles, computing view factors ...")   # progress, before the expensive step
        domain = ViewFactorDomain3D(points, faces, 1, q_in_w, T_in_w, epsilon)   # one element per triangle
        domain(; parallel = true, verbose = false) # analytical view factors between all triangle pairs, multithreaded
        solveEquilibrium!(domain, domain.F_raw; verbose = false)   # single linear solve with the view factors as computed
        println("$T, level $level: largest antipodal asymmetry of F = ",
            Float64(maximum(abs.(domain.F_raw[antipode, antipode] .- domain.F_raw))))   # how far the computed F is from the exact central symmetry

        Tw   = [domain.facesMesh[i].subFaces[1].T_w for i in 1:n_tri]   # surface temperature of every triangle [K]
        free = setdiff(1:n_tri, hot, cold)         # the triangles in radiative equilibrium
        dev  = maximum(abs.(Tw[free] .^ 4 .+ Tw[antipode[free]] .^ 4 .- (T_hot^4 + T_cold^4))) / T_hot^4   # largest relative violation
        println("$T, level $level: largest relative deviation = ", Float64(dev))   # reported in the text
        push!(n_triangles, n_tri)                  # record this level
        push!(deviation, Float64(dev))             # Float64 is enough to plot it (down to ~1e-308)
        finest = domain                            # the last level is the finest
    end
    return n_triangles, deviation, finest
end

## Float64
n_triangles, deviation, finest = symmetry_deviation(Float64, levels)

## 128-bit BigFloat (expensive at level 3, see the cost note in the chat)
setprecision(BigFloat, BITS)                       # working precision of the second series
eps_big = Float64(eps(BigFloat))                   # its machine precision, for the reference line
n_triangles_big, deviation_big, _ = symmetry_deviation(BigFloat, levels_big)

## Figure
fig = Figure(size = (TEXT_WIDTH, 540))             # full text width; height for two panels side by side

scene  = Axis3(fig[1, 1];                          # left: the finest sphere on Makie's standard 3D axis
    xlabel = L"x\,/\,\mathrm{m}",                  # coordinate labels in metres
    ylabel = L"y\,/\,\mathrm{m}",
    zlabel = L"z\,/\,\mathrm{m}",
    aspect = :data,                                # equal scaling, so the sphere looks round
    viewmode = :fitzoom,                           # zoom the 3D box to fill its cell instead of leaving room for every rotation
    azimuth = π / 4,                               # view direction of the LScene default camera (eye at (3, 3, 3)) ...
    elevation = asin(1 / √3),                      # ... looking down at 35°, so the hot cap is visible
    perspectiveness = 0.5,                         # perspective projection, giving depth as in LScene (0 is orthographic, 1 is strongest)
    xypanelvisible = false,                        # no shaded background panels, only the axes and their grid lines,
    yzpanelvisible = false,                        # closer to the open look of LScene
    xzpanelvisible = false,
    xticks = [-1, 0, 1], yticks = [-1, 0, 1], zticks = [-1, 0, 1],   # three ticks per axis suffice for a unit sphere
    xticklabelsize = TICK_SIZE, yticklabelsize = TICK_SIZE, zticklabelsize = TICK_SIZE,   # tick size of the other figures (the theme covers only 2D axes)
    xlabeloffset = 60, ylabeloffset = 60, zlabeloffset = 60)         # move the labels clear of the tick numbers (default 40)
result = plotField(scene, finest; field = :T)      # its surface temperature field (Float64 run), as returned by the package
Colorbar(fig[1, 2], result.plot;                   # the colour bar beside the sphere ...
    label = L"T\,/\,\mathrm{K}")                   # ... labelled with the temperature in kelvin

ax = Axis(fig[1, 3];                               # right: the deviation against the number of triangles
    xscale = log10, yscale = log10,                # both axes logarithmic
    xticks = (n_triangles, string.(n_triangles)),  # tick at every level
    yticks = (10.0 .^ (-40:5:-10), [L"10^{%$k}" for k in -40:5:-10]),   # whole powers of ten, every fifth decade
    limits = (nothing, (1e-40, 1e-10)),            # fixed y-range spanning both machine precisions; x-range automatic
    xlabel = L"\mathrm{Number\ of\ triangles}",
    ylabel = L"\max_k\, |T_k^4 + T_{-k}^4 - T_\mathrm{hot}^4 - T_\mathrm{cold}^4| \,/\, T_\mathrm{hot}^4")
scatterlines!(ax, n_triangles, deviation;          # the largest deviation at each level, Float64
    color = :black, linewidth = 2,                 # black line joining the levels, as in the other figures
    marker = '+', markersize = MARKER_SIZE,        # thin + glyph at the levels
    label = L"\mathrm{Float64}")                   # legend entry
scatterlines!(ax, n_triangles_big, deviation_big;  # the same in 128-bit arithmetic
    color = :black, linewidth = 2,
    marker = '×', markersize = MARKER_SIZE,        # thin × glyph, to tell the series apart without colour
    label = L"\mathrm{128\ bit}")                  # legend entry
hlines!(ax, [eps(Float64), eps_big];               # machine precision of each arithmetic, for reference
    color = :black, linestyle = :dot, linewidth = 2,   # dotted black lines
    label = L"\epsilon_\mathrm{mach}")             # one legend entry for both
Legend(fig[0, 3], ax;                              # legend above the deviation panel, clear of the data as in every figure
    orientation = :horizontal,                     # all entries in one row
    framevisible = false)                          # no frame: it sits outside the axes
colsize!(fig.layout, 1, Aspect(1, 1.0))            # sphere column as wide as row 1 is tall, so no empty space beside the sphere
colgap!(fig.layout, 2, 90)                         # space between the colour bar label and the deviation panel; narrows that panel
fig                                                # show the figure for inspection

## Save, once the figure has been inspected
save("papers/exchange-factor-formulation/figures/sphere_symmetry.pdf", fig)                    # vector output, named after the script
save("papers/exchange-factor-formulation/figures/sphere_symmetry.png", fig; px_per_unit = 2)   # raster copy for inspection, at double resolution