# --- element property text ----------------------------------------------

fmtnum(::Nothing, d, u) = "— " * u
fmtnum(x::Real, d, u)   = string(round(x, digits = d)) * " " * u

# spectral bins are a decomposition of a total — sum them
fmtsum(::Nothing, d, u) = "— " * u
fmtsum(x::Real, d, u)   = string(round(x, digits = d)) * " " * u
fmtsum(x::AbstractVector, d, u) =
    string(round(sum(x), digits = d), " $u (Σ of ", length(x), " bins)")

# per-bin surface property — bins do not add, so report the range
fmtrange(::Nothing, d) = "—"
fmtrange(x::Real, d)   = string(round(x, digits = d))
fmtrange(x::AbstractVector, d) = begin
    lo, hi = extrema(x)
    length(x) == 1 ?
    string(round(x[1], digits = d)) :
    isapprox(lo, hi; atol = 10.0^(-d)) ?
        string(round(lo, digits = d), "  (", length(x), " bins, uniform)") :
        string(round(lo, digits = d), " - ", round(hi, digits = d),
               "  (", length(x), " bins)")
end

faceLoop(sf) = [Point3f(sf.vertices[mod1(j, length(sf.vertices))])
                for j in 1:length(sf.vertices) + 1]

function elementDescription(flat, parent, k)
    sf = flat[k]
    prescribedT = sf.T_in_w > -0.1

    rows = [
        ("area",     "area",                     "$(fmtnum(sf.area, 4,"m²"))"),
        ("midPoint", "midpoint",                 "($(fmtnum(sf.midPoint[1], 3, "m")), $(fmtnum(sf.midPoint[2], 3, "m")), $(fmtnum(sf.midPoint[3], 3, "m")))"),
        ("", "", ""),
        ("epsilon",  "emissivity",               fmtrange(sf.epsilon, 3)),
        ("T_in_w",   "input temperature",        prescribedT ? "$(fmtnum(sf.T_in_w, 1, "K"))" : "unknown"),
        ("q_in_w",   "input source flux power",  prescribedT ? "unknown" : "$(fmtnum(sf.q_in_w, 3, "W"))"),
        ("", "", ""),
        ("T_w",      "temperature",              "$(fmtnum(sf.T_w, 1, "K"))"),
        ("q_w",      "source flux power",        "$(fmtnum(sf.q_w, 3, "W"))"),
        ("j_w",      "radiant power",            "$(fmtsum(sf.j_w, 3, "W"))"),
        ("g_w",      "incident power",           "$(fmtsum(sf.g_w, 3, "W"))"),
        ("g_a_w",    "absorbed incident power",  "$(fmtsum(sf.g_a_w, 3, "W"))"),
        ("e_w",      "emitted radiant power",    "$(fmtsum(sf.e_w, 3, "W"))"),
        ("r_w",      "reflected radiant power",  "$(fmtsum(sf.r_w, 3, "W"))"),
        ("i_w",      "diffuse radiance, j/(πA)", "$(fmtsum(sf.i_w, 3, "W m⁻² sr⁻¹"))"),
    ]

    header = ("FIELD", "DESCRIPTION", "VALUE")
    w1 = maximum(length(r[1]) for r in rows); w1 = max(w1, length(header[1]))
    w2 = maximum(length(r[2]) for r in rows); w2 = max(w2, length(header[2]))
    w3 = maximum(length(r[3]) for r in rows); w3 = max(w3, length(header[3]))
    gap = 2

    line(a, b, c) = rstrip(string(rpad(a, w1 + gap), rpad(b, w2 + gap), c))
    rule = "─"^(w1 + w2 + w3 + 2gap)

    return join(vcat(
        "element $k        superface $(parent[k])",
        "",
        line(header...),
        rule,
        [line(r...) for r in rows],
    ), "\n")
end

# --- unified 3D domain plot ---------------------------------------------

"""
    plotDomain3D(ax, domain; field = :auto, cmap = :thermal, colorrange = nothing,
                 superfaces = false, inspect = false, label = nothing,
                 outlinecolor = :black, linewidth = 4)

Draw `domain` into `ax` as a single mesh whose vertices are ordered by global
element index.

`field = :random` colours each element randomly (mesh inspection); any other
symbol colours by that surface field, e.g. `:T`, `:q`, `:j`. `:auto` picks
`:T` if the domain has been solved and `:random` otherwise. Spectral fields are
summed across bins.

`inspect = true` makes elements hoverable; this requires `DataInspector(fig)`
on the enclosing figure. The hovered element is outlined, and if `label` is a
`Label` placed by the caller, its full property list is written there.

Returns `(plot, selected, describe, elements, parents)`. `plot` is the mesh, so
a colorbar is `Colorbar(fig[1, 2], result.plot)`. `selected` is an
`Observable{Int}` holding the hovered element index.
"""
function plotDomain3D(ax, domain; field = :auto, cmap = :thermal,
                      colorrange = nothing, superfaces::Bool = false,
                      inspect::Bool = false, label = nothing,
                      outlinecolor = :black, linewidth = 4)

    facesVector = domain isa AbstractVector ? domain : domain.facesMesh

    # unmeshed coarse faces: no element indices exist, so no inspection
    if superfaces
        for i in eachindex(facesVector)
            facesVector[i].vertices === nothing && continue
            v = facesVector[i].vertices
            c = RGBf(rand(), rand(), rand())
            Makie.mesh!(ax, [v[1][1], v[2][1], v[3][1]],
                            [v[1][2], v[2][2], v[3][2]],
                            [v[1][3], v[2][3], v[3][3]], color = c)
            if length(v) == 4
                Makie.mesh!(ax, [v[1][1], v[3][1], v[4][1]],
                                [v[1][2], v[3][2], v[4][2]],
                                [v[1][3], v[3][3], v[4][3]], color = c)
            end
        end
        return DomainPlot3D(nothing, nothing, nothing, Any[], Int[])
    end

    flat   = [sf for f in facesVector for sf in f.subFaces]
    parent = [i for (i, f) in enumerate(facesVector) for _ in f.subFaces]

    # choose the colour source
    if field === :auto
        field = flat[1].T_w === nothing ? :random : :T
    end

    if field === :random
        elemcolor = [RGBf(rand(), rand(), rand()) for _ in eachindex(flat)]
        colorkw = (;)
    else
        vals = [getfield(sf, Symbol(string(field) * "_w")) for sf in flat]
        any(isnothing, vals) && error("field :$field is not populated — has the domain been solved?")
        elemcolor = Float32[v isa AbstractVector ? sum(v) : v for v in vals]
        crange = colorrange === nothing ? extrema(elemcolor) : colorrange
        colorkw = (; colormap = cmap, colorrange = crange)
    end

    # one mesh, built in element order; vertices are not shared between
    # elements, so a picked vertex index identifies an element
    nverts = sum(length(sf.vertices) for sf in flat)
    ntris  = sum(length(sf.vertices) == 4 ? 2 : 1 for sf in flat)

    coords    = Matrix{Float32}(undef, nverts, 3)
    conn      = Matrix{Int}(undef, ntris, 3)
    vert2elem = Vector{Int}(undef, nverts)

    vi = 0
    ti = 0
    for (k, sf) in enumerate(flat)
        v = sf.vertices
        base = vi
        for p in v
            vi += 1
            coords[vi, 1] = p[1]
            coords[vi, 2] = p[2]
            coords[vi, 3] = p[3]
            vert2elem[vi] = k
        end
        ti += 1
        conn[ti, 1] = base + 1
        conn[ti, 2] = base + 2
        conn[ti, 3] = base + 3
        if length(v) == 4
            ti += 1
            conn[ti, 1] = base + 1
            conn[ti, 2] = base + 3
            conn[ti, 3] = base + 4
        end
    end

    selected = Makie.Observable(1)
    outline  = Makie.Observable(faceLoop(flat[1]))
    describe(k) = elementDescription(flat, parent, k)

    p = Makie.mesh!(ax, coords, conn; color = elemcolor[vert2elem],
                    inspectable = inspect,
                    inspector_label = (self, i, pos) -> begin
                        1 <= i <= nverts && (selected[] = vert2elem[i])
                        "element $(selected[])"
                    end,
                    colorkw...)

    if inspect
        Makie.lines!(ax, outline; color = outlinecolor, linewidth = linewidth,
                     overdraw = true, inspectable = false)
        Makie.on(selected) do k
            outline[] = faceLoop(flat[k])
            label === nothing || (label.text = describe(k))
        end
        Makie.notify(selected)
    end

    return DomainPlot3D(p, selected, describe, flat, parent)
end

# --- backwards-compatible entry points ----------------------------------

RayTraceHeatTransfer.plotMesh(ax, domain::Union{SurfaceDomain3D{G,P},Vector{PolyFace3D{G}}};
                              field = :random, kwargs...) where {G,P} =
    RayTraceHeatTransfer.plotDomain3D(ax, domain; field = field, kwargs...)

RayTraceHeatTransfer.plotField(ax, domain::Union{SurfaceDomain3D{G,P},Vector{PolyFace3D{G}}};
                               field = :T, kwargs...) where {G,P} =
    RayTraceHeatTransfer.plotDomain3D(ax, domain; field = field, kwargs...)