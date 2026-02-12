# updated version of addSubVolume with spectral support
function addSubVolume!(superVolume::PolyVolume2D{G}, subVolume::PolyVolume2D{G}) where {G}

    # next, update the properties of the added subface
    # the subface inherits the properties of the superface

    # inherit volume state variables - direct assignment/copy
    inheritVolumeProperty!(superVolume, subVolume, :kappa_g)
    inheritVolumeProperty!(superVolume, subVolume, :sigma_s_g)
    inheritVolumeProperty!(superVolume, subVolume, :j_g)
    inheritVolumeProperty!(superVolume, subVolume, :g_a_g)
    inheritVolumeProperty!(superVolume, subVolume, :e_g)
    inheritVolumeProperty!(superVolume, subVolume, :r_g)
    inheritVolumeProperty!(superVolume, subVolume, :g_g)
    inheritVolumeProperty!(superVolume, subVolume, :i_g)
    inheritVolumeProperty!(superVolume, subVolume, :q_in_g)
    inheritVolumeProperty!(superVolume, subVolume, :q_g)
    inheritVolumeProperty!(superVolume, subVolume, :T_in_g)
    inheritVolumeProperty!(superVolume, subVolume, :T_g)

    # then for walls
    for i in 1:length(subVolume.vertices) # loop over walls of subface (same as number of vertices)
        if (length(superVolume.vertices) == 3 && length(subVolume.vertices) == 4) && subVolume.solidWalls[i] == true
            # if superVolume is a triangle and subVolume is a quadrilateral
            # inherits the properties of the superface
            
            # two walls of the quadrilateral need to inherit from the diagonal of the triangle

            # find the diagonal (should still work for equilateral triangles)
            norm1 = norm(superVolume.vertices[1]-superVolume.vertices[2])
            norm2 = norm(superVolume.vertices[2]-superVolume.vertices[3])
            norm3 = norm(superVolume.vertices[3]-superVolume.vertices[1])
            diag_len, diag_index = findmax([norm1, norm2, norm3])
            diag_normal = superVolume.inwardNormals[diag_index]
            
            # check if the outwardNormal of the subface is in the same direction as the diagonal outwardNormal
            diag_dirs = [dot(subVolume.inwardNormals[j], diag_normal) for j in 1:4] # positive for sides which inherit from the diagonal
            # and that the diagonal is solid
            # need to +1 everywhere to account for i=1 is volume
            if diag_dirs[i] > 0.0 && superVolume.solidWalls[diag_index] == true
                # inherits the properties of the superface (the diagonal)
                inheritSurfaceProperties!(superVolume, subVolume; from=diag_index, to=i)

            elseif diag_dirs[i] < 0.0 && superVolume.solidWalls[diag_index] == true
                # if the quadrilateral subface wall is solid, but should not inherit from the superVolume triangle diagonal
                # which index to inherit from? There are four cases
                # 1) subface sides 1 and 2 inherits from diagonal, subface side 3 and 4 inherits from triangle sides 2 and 3, valid
                if diag_dirs[1] > 0.0 && diag_dirs[2] > 0.0
                    # subface side 3 inherits from triangle side 2
                    inheritSurfaceProperties!(superVolume, subVolume; from=2, to=3)

                    # subface side 4 inherits from triangle side 3
                    inheritSurfaceProperties!(superVolume, subVolume; from=3, to=4)
                end

                # 2) subface sides 2 and 3 inherits from diagonal, subface side 1 and 4 inherits from triangle sides 1 and 3, valid
                if diag_dirs[2] > 0.0 && diag_dirs[3] > 0.0
                    # subface side 1 inherits from triangle side 1
                    inheritSurfaceProperties!(superVolume, subVolume; from=1, to=1)

                    # subface side 4 inherits from triangle side 3
                    inheritSurfaceProperties!(superVolume, subVolume; from=3, to=4)
                end

                # 3) subface sides 3 and 4 inherits from diagonal, subface side 1 and 2 inherits from triangle sides 1 and 2, valid
                if diag_dirs[3] > 0.0 && diag_dirs[4] > 0.0
                    # subface side 1 inherits from triangle side 1
                    inheritSurfaceProperties!(superVolume, subVolume; from=1, to=1)

                    # subface side 2 inherits from triangle side 2
                    inheritSurfaceProperties!(superVolume, subVolume; from=2, to=2)                    
                end

                # 4) subface sides 1 and 4 inherits from diagonal, subface side 2 and 3 inherits from triangle sides 1 and 2, valid
                if diag_dirs[4] > 0.0 && diag_dirs[1] > 0.0
                    # subface side 2 inherits from triangle side 1
                    inheritSurfaceProperties!(superVolume, subVolume; from=1, to=2)

                    # subface side 4 inherits from triangle side 3
                    inheritSurfaceProperties!(superVolume, subVolume; from=2, to=3)
                end

            end

        elseif (length(superVolume.vertices) == length(subVolume.vertices)) && subVolume.solidWalls[i] == true
            # if both are quadrilaterals or both are triangles
            inheritSurfaceProperties!(superVolume, subVolume; from=i, to=i)            
        else
            # wall is not solid, no need to inherit
            continue
        end
    end

    # push the subVolume to the list of subVolumes
    push!(superVolume.subVolumes, subVolume)

end