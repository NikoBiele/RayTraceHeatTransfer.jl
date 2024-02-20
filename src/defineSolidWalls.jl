function defineSolidWalls(yLayersHeight::Vector{Float64})
    # This function defines the solid walls of the geometry
    
    solidWalls = Vector{Vector{Bool}}(undef, length(yLayersHeight)-1)

    # below, the walls of each cell are ordered [bottom, right, top, left]
    if length(yLayersHeight) == 2 # then only one sub-enclosure
        solidWalls[1] = [true, true, true, true]
    else
        for i = 1:length(yLayersHeight)-1
            if i == 1
                solidWalls[i] = [true, true, false, true]
            elseif i == length(yLayersHeight)-1
                solidWalls[i] = [false, true, true, true]
            else
                solidWalls[i] = [false, true, false, true]
            end
        end
    end

    return solidWalls
    
end