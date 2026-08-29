# predict_rho_data - deterministic a priori smoothing rate estimate from geometry and F_raw.
#
# Groups surface elements into macro-surfaces (same oriented line), determines
# pairwise blindness by using F_raw, and solves max-weight independent set
# exactly on the small group graph, and returns
#   t_geo   = (2 w(S) - w(V)) / w(V)   (always <= 0 for a valid enclosure,
#                                       by the conservation bound)
#   rho_est = 1 + t_geo/2              (first-order exact near the boundary,
#                                       conservative away from it)
# plus the blind macro-surface set S realizing it.
#
# Gas cells always carry true self-loops (beta > 0) and are assigned to G.
#
# Data-driven a priori AP rate estimate from the traced exchange factors.
# Two macro-surfaces are blind iff no ray connected any pair of their
# elements: the F_raw block between them (either orientation) is all zero.
# This is exactly the coupling the smoother sees, needs no occlusion
# geometry, and has no group-count limit.

function predict_rho_data(rtm::RayTracingDomain2D)
    groups = collect_macro_surfaces(rtm)
    G      = length(groups)
    n_surf = length(rtm.surface_mapping)

    g_of = zeros(Int, n_surf)                      # element id -> group
    for (g, grp) in enumerate(groups), id in grp.element_ids
        g_of[id] = g
    end

    Fbins = rtm.F_raw isa Vector ? rtm.F_raw : [rtm.F_raw]
    wg    = [g.weight for g in groups]

    rho_worst = 0.0; t_worst = -Inf; delta_worst = 0.0
    exact_all = true
    for (b, F) in enumerate(Fbins)
        adj = _surface_adjacency(F, g_of, G, n_surf)

        # a group with internal coupling cannot sit in a blind set
        ok  = [!adj[g, g] for g in 1:G]
        Sg, wS, exact_b = mwis_groups(wg[ok], Matrix(adj[ok, ok]))
        exact_all &= exact_b

        w = get_w(rtm; spectral_bin = b)
        rtm.surfaces_only && (w = w[1:n_surf]; F = F[1:n_surf,1:n_surf])
        W = sum(w)

        t   = (2 * wS - W) / W
        rho_pred = 1 + min(t, 0.0)/2
        if rho_pred > rho_worst
            rho_worst, t_worst = rho_pred, t
            delta_worst = delta_R_raw(F, w)        # measured, not assumed
        end
        wS > W / 2 * (1 + 1e-9) &&
            @warn "predict_rho_data: blind set exceeds half total weight (bin $b)"
    end
    return (; t_data = t_worst, rho_est = rho_worst, delta_init = delta_worst, exact = exact_all)
end

function _surface_adjacency(F::SparseMatrixCSC, g_of, G, n_surf)
    adj  = falses(G, G)
    rows = rowvals(F); vals = nonzeros(F)
    for j in 1:min(size(F, 2), n_surf), k in nzrange(F, j)
        i = rows[k]
        (i <= n_surf && vals[k] > 0) || continue
        a, b = g_of[i], g_of[j]
        adj[a, b] = adj[b, a] = true               # covers both orientations
    end
    return adj
end

function _surface_adjacency(F::Matrix, g_of, G, n_surf)
    adj = falses(G, G)
    for j in 1:n_surf, i in 1:n_surf
        F[i, j] > 0 || continue
        a, b = g_of[i], g_of[j]
        adj[a, b] = adj[b, a] = true
    end
    return adj
end

function collect_macro_surfaces(rtm; tol_n = 1e-9, tol_d = 1e-9)
    raw = []   # (normal, offset, weight, endpoints, id)
    Lref = 0.0
    for ((ci, fi, wi), sidx) in rtm.surface_mapping
        face = rtm.fine_mesh[ci][fi]
        n = face.inwardNormals[wi]
        nn = (n[1], n[2]) ./ hypot(n[1], n[2])
        nverts = length(face.vertices)
        a = face.vertices[wi]
        b = face.vertices[mod1(wi + 1, nverts)]
        d = nn[1]*a[1] + nn[2]*a[2]
        # Orient the normal toward the gas the wall radiates into, using the
        # owning face's midpoint as interior witness. This makes the code
        # independent of the stored normal convention.
        mid = face.midPoint
        if nn[1]*mid[1] + nn[2]*mid[2] - d < 0
            nn = (-nn[1], -nn[2]); d = -d
        end
        push!(raw, (nn, d, face.area[wi],
                    [(a[1], a[2]), (b[1], b[2])], sidx))
        Lref = max(Lref, abs(a[1]), abs(a[2]), abs(b[1]), abs(b[2]))
    end
    Lref = max(Lref, 1.0)
    groups = MacroSurface[]
    for (nn, d, area, eps_, id) in raw
        found = false
        for (g, G) in enumerate(groups)
            if abs(nn[1]*G.normal[1] + nn[2]*G.normal[2] - 1) < tol_n &&
               abs(d - G.offset) < tol_d * Lref
                groups[g] = MacroSurface(G.normal, G.offset, G.weight + area,
                                         vcat(G.endpoints, eps_),
                                         vcat(G.element_ids, id))
                found = true; break
            end
        end
        found || push!(groups, MacroSurface(nn, d, area, eps_, [id]))
    end
    return groups
end

# --- exact max-weight independent set on the group graph ----------
function mwis_groups(wg::Vector{Float64}, adj::Matrix{Bool};
                     node_budget::Int = 1_000_000)
    G = length(wg)
    order = sortperm(wg; rev = true)
    best = Ref(0.0); bestS = Ref(Int[])
    nodes = Ref(0); exact = Ref(true)
    function rec(i, cur, curw, banned)
        if (nodes[] += 1) > node_budget
            exact[] = false
            return
        end
        rem = sum(wg[order[j]] for j in i:G if !banned[order[j]]; init = 0.0)
        curw + rem <= best[] && return
        if i > G
            curw > best[] && (best[] = curw; bestS[] = copy(cur))
            return
        end
        v = order[i]
        if !banned[v]
            b2 = copy(banned)
            for u in 1:G
                adj[v, u] && (b2[u] = true)
            end
            push!(cur, v)
            rec(i + 1, cur, curw + wg[v], b2)
            pop!(cur)
        end
        rec(i + 1, cur, curw, banned)
        curw > best[] && (best[] = curw; bestS[] = copy(cur))
    end
    rec(1, Int[], 0.0, falses(G))
    return bestS[], best[], exact[]
end