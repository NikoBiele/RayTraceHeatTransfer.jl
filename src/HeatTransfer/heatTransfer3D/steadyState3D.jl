function steadyState3D!(domain3D::Domain3D_faces)
    F = domain3D.F
    num_elements = sum([length(domain3D.facesMesh[i].subFaces) for i in 1:length(domain3D.facesMesh)])
    emissivities = zeros(num_elements)
    areas = zeros(num_elements)
    for i in 1:length(domain3D.facesMesh)
        num_subfaces = length(domain3D.facesMesh[i].subFaces)
        for k in 1:num_subfaces
            emissivities[(i-1)*num_subfaces + k] = domain3D.facesMesh[i].subFaces[k].epsilon
            areas[(i-1)*num_subfaces + k] = domain3D.facesMesh[i].subFaces[k].area
        end
    end
    B = ones(length(areas))'*(1.0 .- emissivities)
    K = F .* B
    S_infty = (I-K)\F
    A = (1 .- B)' .* S_infty .* (1 .- B)
    R = (1 .- B)' .* S_infty .* B
    C = I - R' - A'
    D = I - R'
    h = zeros(num_elements) # boundary conditions
    M = zeros(num_elements, num_elements)
    for i in 1:length(domain3D.facesMesh)
        num_subfaces = length(domain3D.facesMesh[i].subFaces)
        for k in 1:num_subfaces
            if domain3D.facesMesh[i].subFaces[k].T_in_w < 0
                h[(i-1)*num_subfaces + k] = domain3D.facesMesh[i].subFaces[k].q_in_w # known flux
                M[(i-1)*num_subfaces + k, :] = C[(i-1)*num_subfaces + k, :]
            else
                h[(i-1)*num_subfaces + k] = domain3D.facesMesh[i].subFaces[k].area*emissivities[(i-1)*num_subfaces + k]*
                                            STEFAN_BOLTZMANN*domain3D.facesMesh[i].subFaces[k].T_in_w^4
                M[(i-1)*num_subfaces + k, :] = D[(i-1)*num_subfaces + k, :]
            end
        end
    end
    j = M \ h # solve the system

    # write results to mesh
    r = R' * j
    abs = A' * j
    g = r + abs
    e = max.(j .- r, 0.0) 
    q = e - abs
    i_w = j ./ (pi .* areas)
    for i in 1:length(domain3D.facesMesh)
        num_subfaces = length(domain3D.facesMesh[i].subFaces)
        for k in 1:num_subfaces
            domain3D.facesMesh[i].subFaces[k].j_w = j[(i-1)*num_subfaces + k]
            domain3D.facesMesh[i].subFaces[k].g_a_w = abs[(i-1)*num_subfaces + k]
            domain3D.facesMesh[i].subFaces[k].e_w = e[(i-1)*num_subfaces + k]
            domain3D.facesMesh[i].subFaces[k].r_w = r[(i-1)*num_subfaces + k]
            domain3D.facesMesh[i].subFaces[k].g_w = g[(i-1)*num_subfaces + k]
            domain3D.facesMesh[i].subFaces[k].i_w = i_w[(i-1)*num_subfaces + k]
            domain3D.facesMesh[i].subFaces[k].q_w = q[(i-1)*num_subfaces + k]
            domain3D.facesMesh[i].subFaces[k].T_w = (e[(i-1)*num_subfaces + k]/
                            (emissivities[(i-1)*num_subfaces + k] * STEFAN_BOLTZMANN * areas[(i-1)*num_subfaces + k]))^(1/4)
        end
    end
    
end