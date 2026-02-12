# used for meshing surfaces only mesh (transparent view factor mesh)
function meshFaces(coords, faces, Ndim)

    num_faces = size(faces, 1)
    points_per_face = size(faces, 2)

    # Project vertices to xy-Plane
    projected_faces = [projectFaceFlat(coords, faces[i, :]) for i in 1:num_faces]

    # mesh each face in the xy-plane
    mesh_pre_project = [points_per_face == 4 ? meshQuad(projected_faces[i][1],Ndim,Ndim) :
                        meshTriangle(projected_faces[i][1],Ndim,Ndim) for i in 1:num_faces]

    # then project each face back
    mesh_post_project = [projectFaceBack(mesh_pre_project[i], projected_faces[i][2], projected_faces[i][3]) for i in 1:num_faces]

    return mesh_post_project
end