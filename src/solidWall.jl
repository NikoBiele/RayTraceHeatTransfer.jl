function solidWall(mesh::TracingMesh, N_subs_count::Int64)
    # this function finds out if any boundaries are not allowed to be crossed
    # if we hit these real walls the ray is absorbed or reflected

    localSolidWalls = mesh.solidWalls[N_subs_count]
    
    return localSolidWalls[1], localSolidWalls[2], localSolidWalls[3], localSolidWalls[4]
end