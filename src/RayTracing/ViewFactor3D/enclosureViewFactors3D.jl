function enclosureViewFactors3D(vfd::ViewFactorDomain3D, parallel::Bool)
    elements = 0
    for superface in vfd.facesMesh
        for subface in superface.subFaces
            elements += 1
        end
    end
            
    F_raw        = zeros(elements, elements)
    faces        = Vector{Matrix{Float64}}(undef, elements)

    # ---- Pass 1: faces vertices --------------------------------
    k = 0
    for superface in vfd.facesMesh
        for subface in superface.subFaces
            k += 1
            verts  = subface.vertices
            face   = Matrix{Float64}(undef, length(verts), 3)
            for v in eachindex(verts)
                face[v, :] = verts[v]
            end
            faces[k] = face
        end
    end

    # ---- Pass 2: upper triangle -------------------------------
    fillRow! = function (i)
        face1 = faces[i]
        for j in (i + 1):elements
            v1, v2, A1, A2 = viewFactor3D(face1, faces[j])
            v1 = (isnan(v1) || v1 < 0.0) ? 0.0 : v1
            v2 = (isnan(v2) || v2 < 0.0) ? 0.0 : v2
            F_raw[i, j] = v1
            F_raw[j, i] = v2
        end
    end

    if parallel
        # Row i costs (elements - i), so a plain @threads split is badly
        # unbalanced. Pairing i with (elements + 1 - i) makes the work per
        # task roughly constant.
        npair = cld(elements, 2)
        @threads for p in 1:npair
            fillRow!(p)
            q = elements + 1 - p
            q > p && fillRow!(q)
        end
    else
        for i in 1:elements
            fillRow!(i)
        end
    end

    return F_raw ./ sum(F_raw, dims=2)
end