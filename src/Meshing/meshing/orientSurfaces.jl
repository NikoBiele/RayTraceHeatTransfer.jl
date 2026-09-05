"""
    faceCycle(row)

Vertex cycle of a face row, dropping the repeated vertex that encodes a
triangle in a 4-wide `faces` matrix. `[1,2,5,5]` → `[1,2,5]`.
"""
function faceCycle(row::AbstractVector{P}) where {P<:Integer}
    v = unique(row)
    length(v) < 3 && error("Degenerate face $(collect(row)): fewer than 3 distinct vertices")
    return v
end

"""
    orientFaces(points, faces) -> (faces_oriented, enclosed_volume)

Rewind every face so the right-hand rule gives the outward normal of the
closed surface. The caller's winding is not trusted.

Faces are first oriented relative to one another by walking the
face-adjacency graph: two faces sharing an edge must traverse it in
opposite directions. That leaves the global sign free, which is then fixed
by requiring the enclosed signed volume to be positive. Neither step
assumes convexity.

Errors if the surface is not a closed manifold — an edge used once is a
hole, an edge used three or more times is a non-manifold junction, and
both would otherwise surface only as escaped rays with no explanation.
"""
function orientFaces(points::AbstractMatrix{G}, faces::AbstractMatrix{P}) where {G,P<:Integer}

    nf     = size(faces, 1)
    cycles = [faceCycle(view(faces, i, :)) for i in 1:nf]

    edgeFaces = Dict{Tuple{P,P},Vector{Int}}()
    for i in 1:nf
        c = cycles[i]
        for k in eachindex(c)
            a = c[k]; b = c[mod1(k + 1, length(c))]
            push!(get!(edgeFaces, a < b ? (a, b) : (b, a), Int[]), i)
        end
    end
    for (key, fs) in edgeFaces
        length(fs) == 2 && continue
        error("Edge $(key) is shared by $(length(fs)) face(s); the enclosure must be a " *
              "closed manifold surface (an edge used once is a hole).")
    end

    hasDirected(c, a, b) =
        any(k -> c[k] == a && c[mod1(k + 1, length(c))] == b, eachindex(c))

    visited = falses(nf); visited[1] = true
    queue   = Int[1]
    while !isempty(queue)
        i  = popfirst!(queue)
        ci = cycles[i]
        for k in eachindex(ci)
            a = ci[k]; b = ci[mod1(k + 1, length(ci))]
            for j in edgeFaces[a < b ? (a, b) : (b, a)]
                (j == i || visited[j]) && continue
                # a consistent neighbour traverses the shared edge as (b, a)
                hasDirected(cycles[j], a, b) && reverse!(cycles[j])
                visited[j] = true
                push!(queue, j)
            end
        end
    end
    all(visited) || error("The surface is disconnected; all faces must form one enclosure")

    o = Point3{G}(vec(sum(points, dims = 1)) ./ size(points, 1))
    V = zero(G)
    for c in cycles
        v1 = Point3{G}(points[c[1], :]) - o
        for k in 2:length(c)-1
            v2 = Point3{G}(points[c[k], :]) - o
            v3 = Point3{G}(points[c[k+1], :]) - o
            V += dot(v1, cross(v2, v3)) / 6
        end
    end
    if V < zero(G)
        for c in cycles
            reverse!(c)
        end
        V = -V
    end

    out  = similar(faces)
    ncol = size(faces, 2)
    for i in 1:nf
        c = cycles[i]
        for k in 1:ncol
            out[i, k] = k <= length(c) ? c[k] : c[end]
        end
    end

    return out, V
end