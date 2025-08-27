function plotField3D(ax1, domain3D::Domain3D_faces, cmap; field=:T)

    # Function to get the field value
    value_field = Symbol(string(field) * "_w")
    
    value_range = extrema([getfield(subface, value_field) for i in 1:length(domain3D.facesMesh) for subface in domain3D.facesMesh[i].subFaces])
    crange = value_range

    num_faces = length(domain3D.facesMesh)
    num_subfaces = length(domain3D.facesMesh[1].subFaces)
    for i in 1:num_faces # face
        for p in 1:num_subfaces
            tri1_x = [domain3D.facesMesh[i].subFaces[p].vertices[1][1], domain3D.facesMesh[i].subFaces[p].vertices[2][1], domain3D.facesMesh[i].subFaces[p].vertices[3][1]]
            tri1_y = [domain3D.facesMesh[i].subFaces[p].vertices[1][2], domain3D.facesMesh[i].subFaces[p].vertices[2][2], domain3D.facesMesh[i].subFaces[p].vertices[3][2]]
            tri1_z = [domain3D.facesMesh[i].subFaces[p].vertices[1][3], domain3D.facesMesh[i].subFaces[p].vertices[2][3], domain3D.facesMesh[i].subFaces[p].vertices[3][3]]
            plot = Makie.mesh!(ax1, tri1_x, tri1_y, tri1_z, 
                                color = getfield(domain3D.facesMesh[i].subFaces[p], value_field),
                                colormap = cmap, 
                                colorrange = crange)
            if length(domain3D.facesMesh[i].subFaces[p].vertices) == 4
                tri2_x = [domain3D.facesMesh[i].subFaces[p].vertices[1][1], domain3D.facesMesh[i].subFaces[p].vertices[3][1], domain3D.facesMesh[i].subFaces[p].vertices[4][1]]
                tri2_y = [domain3D.facesMesh[i].subFaces[p].vertices[1][2], domain3D.facesMesh[i].subFaces[p].vertices[3][2], domain3D.facesMesh[i].subFaces[p].vertices[4][2]]
                tri2_z = [domain3D.facesMesh[i].subFaces[p].vertices[1][3], domain3D.facesMesh[i].subFaces[p].vertices[3][3], domain3D.facesMesh[i].subFaces[p].vertices[4][3]]
                plot = Makie.mesh!(ax1, tri2_x, tri2_y, tri2_z, 
                            color = getfield(domain3D.facesMesh[i].subFaces[p], value_field),
                            colormap = cmap, 
                            colorrange = crange)
            end
        end
    end

end