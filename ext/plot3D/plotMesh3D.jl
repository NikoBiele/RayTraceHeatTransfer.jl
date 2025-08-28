function RayTraceHeatTransfer.plotMesh3D(ax1, domain3D::Domain3D_faces)

    # Determine the "grid" dimensions - assuming all faces have same structure
    num_faces = length(domain3D.facesMesh)
    num_subFaces = length(domain3D.facesMesh[1].subFaces)  # Total points in first sub-array
    # For your data, this appears to be 15 points arranged in a triangular pattern
    
    facecolor = rand(num_faces, num_subFaces, 3)
    color_save = zeros(num_faces, num_subFaces, 3)

    for i in 1:num_faces # face
        for k in 1:num_subFaces # subface

            color_save[i, k, :] .= facecolor[i, k, :]
            
            # Create triangles using the quad logic, but with linear indexing
            tri1_x = [domain3D.facesMesh[i].subFaces[k].vertices[1][1], domain3D.facesMesh[i].subFaces[k].vertices[2][1], domain3D.facesMesh[i].subFaces[k].vertices[3][1]]
            tri1_y = [domain3D.facesMesh[i].subFaces[k].vertices[1][2], domain3D.facesMesh[i].subFaces[k].vertices[2][2], domain3D.facesMesh[i].subFaces[k].vertices[3][2]]
            tri1_z = [domain3D.facesMesh[i].subFaces[k].vertices[1][3], domain3D.facesMesh[i].subFaces[k].vertices[2][3], domain3D.facesMesh[i].subFaces[k].vertices[3][3]]
            Makie.mesh!(ax1, tri1_x, tri1_y, tri1_z, color = RGBf(color_save[i, k, :]...))
            if length(domain3D.facesMesh[i].subFaces[k].vertices) == 4
                tri2_x = [domain3D.facesMesh[i].subFaces[k].vertices[1][1], domain3D.facesMesh[i].subFaces[k].vertices[3][1], domain3D.facesMesh[i].subFaces[k].vertices[4][1]]
                tri2_y = [domain3D.facesMesh[i].subFaces[k].vertices[1][2], domain3D.facesMesh[i].subFaces[k].vertices[3][2], domain3D.facesMesh[i].subFaces[k].vertices[4][2]]
                tri2_z = [domain3D.facesMesh[i].subFaces[k].vertices[1][3], domain3D.facesMesh[i].subFaces[k].vertices[3][3], domain3D.facesMesh[i].subFaces[k].vertices[4][3]]
                Makie.mesh!(ax1, tri2_x, tri2_y, tri2_z, color = RGBf(color_save[i, k, :]...))
            end
        end
    end
end