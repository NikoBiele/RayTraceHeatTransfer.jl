function RayTracingDomain3D_surfaces(points::Matrix{G}, faces::Matrix{P}, Ndims::P,
                                     q_in_w::Vector{G}, T_in_w::Vector{G},
                                     epsilon::Union{Vector{G}, Vector{Vector{G}}};
                                     planarity_tol = 1e-8) where {G, P<:Integer}

    is_spectral = isa(epsilon[1], Vector)
    if is_spectral
        if std([std(epsilon[i]) for i in 1:size(epsilon, 1)]) > 1e-6
            spectral_mode = :spectral_variable
        else
            spectral_mode = :spectral_uniform
        end
    else
        spectral_mode = :grey
    end
    n_bins = is_spectral ? length(epsilon[1]) : 1

    faces_oriented, enclosed_volume = orientFaces(points, faces)

    # domain_midpoint = nothing → normals from winding, valid for non-convex
    superFaces, uniform_epsilon =
        buildSurfaceMesh(points, faces_oriented, Ndims, epsilon, q_in_w, T_in_w, n_bins;
                         planarity_tol = planarity_tol)

    return RayTracingDomain3D_surfaces{G,P}(
        points, faces_oriented, Ndims, superFaces,
        spzeros(G, 0, 0), spzeros(G, 0, 0),
        nothing, nothing, nothing, false,
        spectral_mode, n_bins, nothing, nothing, uniform_epsilon, true)
end

function buildSurfaceMesh(points::AbstractMatrix, faces::AbstractMatrix, Ndims::Integer,
                          epsilon, q_in_w, T_in_w, n_bins::Int;
                          planarity_tol = 1e-8)

    G = eltype(points)

    superFaces = PolyFace3D[]
    for (i, face_rows) in enumerate(eachrow(faces))
        if length(unique(face_rows)) == 3
            points3d = [Point3{G}(points[fr, :]) for fr in unique(face_rows)]
        elseif length(face_rows) == 4
            points3d = [Point3{G}(points[face_rows[j], :]) for j in 1:4]
        else
            error("Superfaces in 3D must consist of 3 or 4 points")
        end
        checkPlanarFace(points3d, i, G(planarity_tol))
        inwardNormal = -normalize(cross(points3d[2] - points3d[1], points3d[3] - points3d[1]))
        push!(superFaces, PolyFace3D(points3d, true, inwardNormal,
                                     epsilon[i], q_in_w[i], T_in_w[i]))
    end

    mesh3D = meshFaces(points, faces, Ndims)

    first_epsilon = nothing
    uniform_epsilon = true
    for i in 1:length(mesh3D)
        superFaces[i].subFaces = PolyFace3D[]
        for j in 1:length(mesh3D[i][1])
            p1 = Point3{G}(mesh3D[i][1][j]...)
            p2 = Point3{G}(mesh3D[i][2][j]...)
            p3 = Point3{G}(mesh3D[i][3][j]...)
            p4 = Point3{G}(mesh3D[i][4][j]...)
            inwardNormal = -normalize(cross(p2 - p1, p3 - p1))
            vlist = isapprox(p3, p4, atol = 1e-5) ? [p1, p2, p3] : [p1, p2, p3, p4]
            push!(superFaces[i].subFaces,
                  PolyFace3D(vlist, true, inwardNormal, epsilon[i], 0.0, T_in_w[i]))
        end

        total_area = sum(sf.area for sf in superFaces[i].subFaces)
        for subface in superFaces[i].subFaces
            # Each subface gets flux proportional to its area fraction
            subface.q_in_w = q_in_w[i] * (subface.area / total_area)
            if first_epsilon === nothing && isa(subface.epsilon, Vector)
                first_epsilon = subface.epsilon[1]
            elseif first_epsilon === nothing && !isa(subface.epsilon, Vector)
                first_epsilon = subface.epsilon
            else
                for bin in 1:n_bins
                    if !isapprox(subface.epsilon[bin], first_epsilon, atol=1e-5)
                        uniform_epsilon = false
                    end
                end
            end
        end
    end

    return superFaces, uniform_epsilon
end

function (dom::RayTracingDomain3D_surfaces{G,P})(rays_tot::Integer; nudge = nothing,
                            nthreads::S=Threads.nthreads(),
                            seeds::Union{Vector{K},K,UnitRange{K}}=1:Threads.nthreads(),
                            rngs::Union{Vector{<:Random.AbstractRNG},<:Random.AbstractRNG}=
                                Random.Xoshiro.(1:Threads.nthreads()),
                            verbose::Bool = true) where {G,P<:Integer,K<:Integer,S<:Integer}
    trace_nudge = nudge === nothing ? G(10_000) * eps(G) : G(nudge)
    verbose && println("Ray tracing surface enclosure " *
                       "(geometry only, wavelength-independent)...")
    traceSurfaces3D(dom, rays_tot, trace_nudge; nthreads=nthreads,
                    seeds=seeds, rngs=rngs, verbose=verbose)
    return nothing
end