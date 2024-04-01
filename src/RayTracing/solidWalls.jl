"""
    solidWalls(mesh::RayTracingMesh, N_subs_count::Int64)

This function finds out if any boundaries are not allowed to be crossed.
If we hit these real walls the ray is absorbed or reflected.
"""
function solidWalls(mesh::RayTracingMesh, N_subs_count::Int64)
    
    localSolidWalls1 = mesh.solidWalls[N_subs_count][1]
    localSolidWalls2 = mesh.solidWalls[N_subs_count][2]
    localSolidWalls3 = mesh.solidWalls[N_subs_count][3]
    localSolidWalls4 = mesh.solidWalls[N_subs_count][4]

    return localSolidWalls1, localSolidWalls2, localSolidWalls3, localSolidWalls4
end