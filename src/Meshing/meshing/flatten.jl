function orthonormalFrame(n::Point3{G}) where {G}
    a  = abs(n[1]) < 0.9 ? Point3{G}(1, 0, 0) : Point3{G}(0, 1, 0)
    e1 = normalize(cross(a, n))
    e2 = cross(n, e1)
    return e1, e2
end

function flatten!(domain::RayTracingDomain3D_surfaces{G}) where {G}
    facets = Facet3D{G}[]
    tris   = Tri3D{G}[]

    for (p, superFace) in enumerate(domain.facesMesh)
        leaves = (superFace.subFaces === nothing || isempty(superFace.subFaces)) ?
                 [superFace] : superFace.subFaces

        for leaf in leaves
            idx = Int32(length(facets) + 1)
            v   = leaf.vertices
            inwardNormal = -normalize(cross(v[2] - v[1], v[3] - v[1]))
            nv  = Int8(length(v))
            e1, e2 = orthonormalFrame(inwardNormal)

            if nv == 3
                vv    = SVector(v[1], v[2], v[3], v[3])
                split = one(G)
            else
                vv = SVector(v[1], v[2], v[3], v[4])
                a1 = norm(cross(v[2] - v[1], v[3] - v[1])) / 2
                a2 = norm(cross(v[3] - v[1], v[4] - v[1])) / 2
                split = a1 / (a1 + a2)
            end

            push!(facets, Facet3D{G}(vv, nv, inwardNormal, e1, e2,
                                     Point3{G}(leaf.midPoint), leaf.area, split, Int32(p)))

            push!(tris, Tri3D{G}(vv[1], vv[2] - vv[1], vv[3] - vv[1], idx))
            nv == 4 && push!(tris, Tri3D{G}(vv[1], vv[3] - vv[1], vv[4] - vv[1], idx))
        end
    end

    domain.facets    = facets
    domain.tris      = tris
    domain.bvh       = buildBVH(tris)
    domain.flattened = true
    return domain
end