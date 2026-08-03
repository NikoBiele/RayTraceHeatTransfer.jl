# Updated constructor from existing RayTracingMesh - now with spectral support
function RayTracingDomain2D(rtm::IntermediateMesh2D, verbose::Bool)
    # Extract the type parameters from the input
    VPF = typeof(rtm.coarse_mesh)
    VVPF = typeof(rtm.fine_mesh)
    MT = rtm.F_raw isa Vector ? typeof(rtm.F_raw[1]) : typeof(rtm.F_raw)
    GRID = typeof(rtm.coarse_grid)
    
    # Build optimized cache structures
    
    # 1. Coarse face cache (direct vector access)
    coarse_face_cache = [face for face in rtm.coarse_mesh]
    
    # 2. Fine face cache (flattened structure) 
    fine_face_cache = [Vector{eltype(submesh)}([face for face in submesh]) for submesh in rtm.fine_mesh]
    
    # 3. Pre-compute geometric data for coarse mesh
    coarse_wall_normals = [copy(face.inwardNormals) for face in rtm.coarse_mesh]
    coarse_wall_midpoints = [copy(face.wallMidPoints) for face in rtm.coarse_mesh]
    
    coarse_bounding_boxes = Vector{Tuple{Point2, Point2}}(undef, length(rtm.coarse_mesh))
    for (i, face) in enumerate(rtm.coarse_mesh)
        vertices = face.vertices
        min_x = minimum(v[1] for v in vertices)
        max_x = maximum(v[1] for v in vertices)
        min_y = minimum(v[2] for v in vertices)
        max_y = maximum(v[2] for v in vertices)
        coarse_bounding_boxes[i] = (Point2(min_x, min_y), Point2(max_x, max_y))
    end
    
    # 4. Pre-compute geometric data for fine mesh
    fine_wall_normals = Vector{Vector{Vector{Point2}}}(undef, length(rtm.fine_mesh))
    fine_wall_midpoints = Vector{Vector{Vector{Point2}}}(undef, length(rtm.fine_mesh))
    fine_bounding_boxes = Vector{Vector{Tuple{Point2, Point2}}}(undef, length(rtm.fine_mesh))
    
    for (i, submesh) in enumerate(rtm.fine_mesh)
        fine_wall_normals[i] = [copy(face.inwardNormals) for face in submesh]
        fine_wall_midpoints[i] = [copy(face.wallMidPoints) for face in submesh]
        
        fine_bounding_boxes[i] = Vector{Tuple{Point2, Point2}}(undef, length(submesh))
        for (j, face) in enumerate(submesh)
            vertices = face.vertices
            min_x = minimum(v[1] for v in vertices)
            max_x = maximum(v[1] for v in vertices)
            min_y = minimum(v[2] for v in vertices)
            max_y = maximum(v[2] for v in vertices)
            fine_bounding_boxes[i][j] = (Point2(min_x, min_y), Point2(max_x, max_y))
        end
    end
    
    # 5. Determine spectral mode from the first face
    first_face = fine_face_cache[1][1]
    is_spectral = isa(first_face.kappa_g, Vector)
    n_spectral_bins = is_spectral ? length(first_face.kappa_g) : 1
    
    # 6. Compute mappings
    surface_mapping = Dict{Tuple{Int,Int,Int}, Int}()
    volume_mapping = Dict{Tuple{Int,Int}, Int}()
    surface_index = 1
    volume_index = 1
    surface_areas = []
    volumes = []
    for (coarse_index, coarse_face) in enumerate(rtm.coarse_mesh)
        for (fine_index, fine_face) in enumerate(rtm.fine_mesh[coarse_index])
            for (wall_index, is_solid) in enumerate(fine_face.solidWalls)
                if is_solid
                    surface_mapping[(coarse_index, fine_index, wall_index)] = surface_index
                    push!(surface_areas, fine_face.area[wall_index])
                    surface_index += 1
                end
            end
            volume_mapping[(coarse_index, fine_index)] = volume_index
            push!(volumes, fine_face.volume)
            volume_index += 1
        end
    end
    VT = typeof(surface_areas)
    DIII = typeof(surface_mapping)
    DII = typeof(volume_mapping)
        
    rtm_optim = RayTracingDomain2D{VPF,VVPF,MT,VT,DIII,DII,GRID}(
        rtm.coarse_mesh, rtm.fine_mesh, rtm.coarse_grid, rtm.fine_grids,
        rtm.F_raw, rtm.F_smooth,
        surface_areas, volumes, surface_mapping, volume_mapping,
        :not_yet_set, n_spectral_bins, nothing, # NEW spectral fields
        false, [0.0],
        coarse_face_cache, fine_face_cache,
        coarse_wall_normals, coarse_wall_midpoints,
        fine_wall_normals, fine_wall_midpoints,
        coarse_bounding_boxes, fine_bounding_boxes,
        # Initialize spatial acceleration as nothing (build later)
        nothing,  # coarse_grid_opt
        nothing,  # coarse_bboxes_opt  
        nothing,  # fine_grids_opt
        nothing,  # fine_bboxes_opt
        nothing   # energy_error
    )

    rtm_optim.uniform_across_bin = validateExtinctionUniformity!(rtm_optim; verbose=verbose)

    # Update spectral mode based on uniformity
    if validateSpectralUniformity!(rtm_optim; verbose=verbose) && is_spectral
        rtm_optim.spectral_mode = :spectral_uniform
    elseif is_spectral
        rtm_optim.spectral_mode = :spectral_variable
    else
        rtm_optim.spectral_mode = :grey
    end

    return rtm_optim
end

# Updated constructor that builds from scratch (for new meshes) - now with spectral support  
function RayTracingDomain2D(faces::Vector{PolyVolume2D{G}}, Ndiv::Vector{Tuple{P,P}};
                            verbose::Bool=true) where {G, P<:Integer}
    # First create the standard RayTracingMesh
    verbose && println("Building intermediate mesh...")
    standardMesh = IntermediateMesh2D(faces, Ndiv)
    
    # Then convert to optimized version with spectral support
    verbose && println("Optimizing mesh...")
    optimMesh = RayTracingDomain2D(standardMesh, verbose)

    surfaces_only = true
    for face in faces
        if face.volume*sum(face.kappa_g + face.sigma_s_g)/length(face.kappa_g) > 1e-8
            surfaces_only = false
            break
        end
    end
    optimMesh.surfaces_only = surfaces_only

    # Build spatial acceleration
    verbose && println("Building spatial acceleration structures...")
    buildSpatialAcceleration!(optimMesh)

    # predict smoothing rate based on geometry (detect parallel plate geometries)
    predict_geom = predict_rho_geometric(optimMesh)

    # warn the user if parallel plates are detected
    N = length(optimMesh.surface_mapping) + length(optimMesh.volume_mapping)
    guard = sqrt(N) * 8 * eps(Float64)
    delta_init_est = 0.1 # estimate of initial violation
    iters_required = trunc(Int, log(guard/delta_init_est)/log(predict_geom.rho_est))
    if predict_geom.rho_est > 0.99
        @warn "Parallel-plate-like geometry detected.\n"*
        "Convergence rate for pure alternating-projection (AP) smoothing predicted to be ρ ≈ $(round(predict_geom.rho_est; digits=8)).\n"*
        "It is recommended to use pure Dykstra iteration -if feasible- (it produces dense matrices, and is slow at scale).\n"*
        "Call for example like: 'mesh(N_rays; method=:exchange, k_dykstra=1000, max_iters=0)'.\n"*
        "This enables 1000 Dykstra rounds and disables AP. Alternatively, trace more rays and use 'F_raw'.\n"*
        "If Dykstra is infeasible, increase AP 'max_iters' to at least $iters_required for succesfull AP smoothing.\n"
    end

    return optimMesh
end

# predict_rho_geometric -- deterministic a priori rate estimate from geometry.
#
#   report = predict_rho_geometric(rtm)          # rtm::RayTracingDomain2D
#
# Groups surface elements into macro-surfaces (same oriented line), determines
# pairwise blindness by exact geometric tests (collinear / facing-away /
# occluded), solves max-weight independent set EXACTLY on the small group
# graph, and returns
#   t_geo   = (2 w(S) - w(V)) / w(V)   (always <= 0 for a valid enclosure,
#                                       by the conservation bound)
#   rho_est = 1 + t_geo/2              (first-order exact near the boundary,
#                                       conservative away from it)
# plus the blind macro-surface set S realizing it.
#
# Gas cells always carry true self-loops (beta > 0) and are assigned to G.
# Complexity: O(N_surf) grouping + O(G^2 * N_solid) blindness + exact MWIS on
# G groups (G = number of distinct wall lines; refinement-invariant).

#----------- rate and t prediction from pure geometry --------------

# strict segment intersection (shared endpoints do not count as blocking)
function segs_intersect(p1, p2, q1, q2; eps = 1e-12)
    d1 = p2 - p1; d2 = q2 - q1
    den = d1[1]*d2[2] - d1[2]*d2[1]
    abs(den) < eps && return false
    dq = q1 - p1
    s = (dq[1]*d2[2] - dq[2]*d2[1]) / den
    t = (dq[1]*d1[2] - dq[2]*d1[1]) / den
    return (eps < s < 1 - eps) && (eps < t < 1 - eps)
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

# ---------- pairwise blindness between macro-surfaces ----------

function blind_pair(A::MacroSurface, B::MacroSurface, solid_segs, G; tol = 1e-9)
    # (1) same oriented line -> collinear, exactly blind
    if abs(A.normal[1]*B.normal[1] + A.normal[2]*B.normal[2] - 1) < tol &&
       abs(A.offset - B.offset) < tol
        return true
    end
    # (2) facing away: B entirely in A's back half-space, or vice versa
    behindA = all(p -> A.normal[1]*p[1] + A.normal[2]*p[2] - A.offset <= tol,
                  B.endpoints)
    behindB = all(p -> B.normal[1]*p[1] + B.normal[2]*p[2] - B.offset <= tol,
                  A.endpoints)
    (behindA || behindB) && return true
    # (3) occlusion: every connecting segment blocked by some solid wall
    if G <= 100
        own = Set(vcat(A.element_ids, B.element_ids))
        for pa in A.endpoints, pb in B.endpoints
            blocked = false
            for (q1, q2, id) in solid_segs
                id in own && continue
                if segs_intersect([pa[1], pa[2]], [pb[1], pb[2]],
                                [q1[1], q1[2]], [q2[1], q2[2]])
                    blocked = true; break
                end
            end
            blocked || return false      # one unblocked sight line -> visible
        end
    end
    return true
end

# ---------- exact max-weight independent set on the group graph ----------

function mwis_groups(wg::Vector{Float64}, adj::Matrix{Bool})
    G = length(wg)
    order = sortperm(wg; rev = true)
    best = Ref(0.0); bestS = Ref(Int[])
    function rec(i, cur, curw, banned)
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
    return bestS[], best[]
end

# ---------- main entry ----------

function predict_rho_geometric(rtm)
    groups = collect_macro_surfaces(rtm)
    G = length(groups)
    if G > 100
        @warn "Too many macro surfaces in the domain to include the occlusion test: "*
              "Domain may or may not be parallel-plate-like, and alternating-projection smoothing might stall."
    end

    # solid wall segments for occlusion tests, tagged by element id
    solid_segs = []
    for g in groups, k in 1:2:length(g.endpoints)
        push!(solid_segs, (g.endpoints[k], g.endpoints[k+1],
                           g.element_ids[div(k + 1, 2)]))
    end

    adj = falses(G, G)               # edge = potentially visible (NOT blind)
    for a in 1:G, b in (a+1):G
        adj[a, b] = adj[b, a] = !blind_pair(groups[a], groups[b], solid_segs, G)
    end

    wg = [g.weight for g in groups]
    Sg, wS = mwis_groups(wg, Matrix(adj))

    w = RayTraceHeatTransfer.get_w(rtm)
    W = sum(w)                        # gas included in the total, always in G
    t = (2 * wS - W) / W

    # Conservation guard: true geometry cannot contain a blind set past half
    # weight. Tripping this means a blindness test mis-fired (orientation,
    # tolerance, or occlusion edge case) -- distrust the result and inspect.
    if wS > W / 2 * (1 + 1e-9)
        @warn "predict_rho_geometric: blind set exceeds half total weight " *
              "(wS=$wS, W=$W). This violates the conservation bound for " *
              "valid enclosures; a geometric blindness test likely mis-fired."
    end

    return (; t_geo = t,
            rho_est = 1 + min(t, 0.0) / 2)
end