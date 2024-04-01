"""
    defineSolidWalls(subs::Vector{SubEnclosure})

This function defines the solid walls of the mesh based on the user inputs for each SubEnclosure.
"""
function defineSolidWalls(subs::Vector{SubEnclosure})
    
    solidWalls = Vector{Vector{Bool}}(undef, length(subs))

    # below, the walls of each cell are ordered [bottom, right, top, left]
    if length(subs) == 1 # then only one sub-enclosure
        solidWalls[1] = Vector{Bool}(undef, 4)
        solidWalls[1] .= [true, true, true, true]
    else
        for i = 1:length(subs)
            solidWalls[i] = Vector{Bool}(undef, 4)
            solidWalls[i] .= [subs[i].solidWall1, subs[i].solidWall2, subs[i].solidWall3, subs[i].solidWall4]
        end
    end

    return solidWalls
    
end