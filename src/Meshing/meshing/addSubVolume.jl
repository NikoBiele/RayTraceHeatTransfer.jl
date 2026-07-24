# updated version of addSubVolume with spectral support
function addSubVolume!(superVolume::PolyVolume2D{G}, subVolume::PolyVolume2D{G}) where {G}

    # next, update the properties of the added subface
    # the subface inherits the properties of the superface

    # inherit volume state variables - direct assignment/copy
    inheritVolumeProperty!(superVolume, subVolume, :kappa_g)
    inheritVolumeProperty!(superVolume, subVolume, :sigma_s_g)
    inheritVolumeProperty!(superVolume, subVolume, :j_g)
    inheritVolumeProperty!(superVolume, subVolume, :g_a_g)
    inheritVolumeProperty!(superVolume, subVolume, :e_g)
    inheritVolumeProperty!(superVolume, subVolume, :r_g)
    inheritVolumeProperty!(superVolume, subVolume, :g_g)
    inheritVolumeProperty!(superVolume, subVolume, :i_g)
    inheritVolumeProperty!(superVolume, subVolume, :q_in_g)
    inheritVolumeProperty!(superVolume, subVolume, :q_g)
    inheritVolumeProperty!(superVolume, subVolume, :T_in_g)
    inheritVolumeProperty!(superVolume, subVolume, :T_g)

    # walls: a solid subface wall lies on exactly one superface edge; find it
    # geometrically and inherit from that edge. Works for tri->quad, tri->tri,
    # quad->quad, any aspect ratio, no diagonal heuristic.
    charlen = maximum(norm(superVolume.vertices[k] - superVolume.vertices[mod1(k+1, length(superVolume.vertices))])
                      for k in 1:length(superVolume.vertices))
    for i in 1:length(subVolume.vertices)
        subVolume.solidWalls[i] || continue
        m = (subVolume.vertices[i] + subVolume.vertices[mod1(i + 1, length(subVolume.vertices))]) / 2
        k, d = containing_edge(superVolume, m)
        if d < 1e-8 * charlen && superVolume.solidWalls[k]
            inheritSurfaceProperties!(superVolume, subVolume; from = k, to = i)
        elseif d >= 1e-8 * charlen
            @warn "Solid subface wall midpoint not on any superface edge (d=$d); no inheritance" 
        end
    end

    # push the subVolume to the list of subVolumes
    push!(superVolume.subVolumes, subVolume)

end

# which superface edge does point p lie on? (edge k spans vertices[k] -> vertices[k+1])
function containing_edge(superVolume::PolyVolume2D, p)
    nV = length(superVolume.vertices)
    best_d, best_k = Inf, 0
    for k in 1:nV
        a = superVolume.vertices[k]
        b = superVolume.vertices[mod1(k + 1, nV)]
        ab = b - a
        t  = clamp(dot(p - a, ab) / dot(ab, ab), 0.0, 1.0)
        d  = norm(p - (a + t * ab))
        if d < best_d
            best_d, best_k = d, k
        end
    end
    return best_k, best_d
end