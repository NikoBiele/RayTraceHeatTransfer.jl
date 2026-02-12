# Updated constructor
function IntermediateMesh2D(faces::Vector{PolyVolume2D{G}}, Ndiv::Vector{Tuple{P,P}}) where {G, P<:Integer}
    meshing_mesh = PolyVolume2D{G}[]
    for (index, face) in enumerate(faces)
        if length(face.vertices) == 3
            if Ndiv[index][1] != Ndiv[index][2]
                error("Number of divisions must be equal for triangles.")
            else
                push!(meshing_mesh, meshTriangle(face, Ndiv[index][1]))
            end
        elseif length(face.vertices) == 4
            push!(meshing_mesh, meshQuad(face, Ndiv[index][1], Ndiv[index][2]))
        else
            error("Only triangles and quadrilaterals are supported.")
        end
    end
    
    coarse_mesh = faces # unmeshed
    if length(meshing_mesh) == 1
        fine_mesh = [meshing_mesh[1].subVolumes] # meshed
    else
        fine_mesh = [submesh.subVolumes for submesh in meshing_mesh] # meshed
    end
    
    # Build grids for coarse and fine meshes
    coarse_grid = buildSpatialGrid(coarse_mesh, one(P))
    fine_grids = [buildSpatialGrid(submesh, one(P)) for submesh in fine_mesh]
    
    # Determine spectral mode from the first face
    first_face = fine_mesh[1][1]
    is_spectral = isa(first_face.kappa_g, Vector)
    
    # Initialize F matrices based on spectral mode
    if is_spectral
        # Spectral mode - create vector of matrices (will be populated during ray tracing)
        F_raw = Matrix{G}[]
        push!(F_raw, zeros(G, 2, 2))
        F_smooth = Matrix{G}[]
        push!(F_smooth, zeros(G, 2, 2))
    else
        # Grey mode - single matrices
        F_raw = zeros(G, 2, 2)
        F_smooth = zeros(G, 2, 2)
    end
    
    rtm = IntermediateMesh2D(
        coarse_mesh, 
        fine_mesh,
        coarse_grid,
        fine_grids,
        F_raw,
        F_smooth,
    )

    return rtm
end