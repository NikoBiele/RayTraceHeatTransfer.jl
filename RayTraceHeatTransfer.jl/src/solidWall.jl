function solidWall(Nx::Int64, Ny::Int64, xCount::Int64, yCount::Int64, N_subs::Int64)
    # this function finds out if any boundaries are not allowed to be crossed
    # if we hit these real walls the ray is absorbed or reflected
    # therefore we set a counter if this is the case
    
    solidWallx = 0
    if xCount == 1
        # left wall is the outer boundary
        solidWallx = 4
    elseif xCount == Nx
        # right wall is the outer boundary
        solidWallx = 2
    end
    solidWally = 0
    if yCount == 1
        # bottom wall is the outer boundary
        solidWally = 1
    elseif yCount == Ny*N_subs
        # top wall is the outer boundary
        solidWally = 3
    end
    
    return solidWallx, solidWally
end