# used for meshing surfaces only mesh (transparent view factor mesh)
function meshFaces(coords, faces, Ndim)

    num_faces = size(faces, 1)

    # Project vertices to xy-Plane
    projected_faces = []
    for i in 1:num_faces
        face_rows = faces[i,:]
        if length(unique(face_rows)) == 3 || length(face_rows) == 3
            push!(projected_faces, projectFaceFlat(coords, unique(face_rows)))
        elseif length(face_rows) == 4
            push!(projected_faces, projectFaceFlat(coords, faces[i, 1:4]))
        end
    end    

    # mesh each face in the xy-plane
    mesh_pre_project = []
    for i in 1:num_faces
        face_rows = faces[i,:]
        if length(unique(face_rows)) == 3 || (length(face_rows) == 3)
            push!(mesh_pre_project, meshTriangle(projected_faces[i][1],Ndim))
        elseif length(face_rows) == 4
            push!(mesh_pre_project, meshQuad(projected_faces[i][1],Ndim,Ndim))
        end
    end

    # then project each face back
    mesh_post_project = [projectFaceBack(mesh_pre_project[i], projected_faces[i][2], projected_faces[i][3]) for i in 1:num_faces]

    return mesh_post_project
end