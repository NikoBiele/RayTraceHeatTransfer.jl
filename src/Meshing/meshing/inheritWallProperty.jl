# Helper function to inherit wall properties
function inheritWallProperty!(superFace::PolyVolume2D{G}, subFace::PolyVolume2D{G}, 
                               property::Symbol, from::Int, to::Int) where {G}
    if property == :q_in_w || property == :q_w
        
        super_wall_array = getfield(superFace, property)
        sub_wall_array = getfield(subFace, property)
        
        super_val = super_wall_array[from]
        sub_val = sub_wall_array[to]
        
        # source flux should reduce by area ratio
        sub_wall_array[to] = super_val*subFace.area[to]/superFace.area[from]
    else
        super_wall_array = getfield(superFace, property)
        sub_wall_array = getfield(subFace, property)
        
        super_val = super_wall_array[from]
        sub_val = sub_wall_array[to]
        
        if isa(super_val, Vector)
            # Spectral - copy entire vector
            sub_val .= super_val
        else
            # Grey - direct assignment
            sub_wall_array[to] = super_val
        end
    end
end