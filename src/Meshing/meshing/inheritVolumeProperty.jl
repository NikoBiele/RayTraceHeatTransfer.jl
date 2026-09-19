# Helper function to inherit volume properties
function inheritVolumeProperty!(superFace::PolyVolume2D{G}, subFace::PolyVolume2D{G}, property::Symbol) where {G}
    if property == :phase
        # descriptors are immutable: assign directly (a per-bin vector is copied)
        super_val = getfield(superFace, property)
        setfield!(subFace, property, super_val isa Vector ? copy(super_val) : super_val)
    elseif property == :q_in_g || property == :q_g
        super_val = getfield(superFace, property)
        sub_val = getfield(subFace, property)

        # source flux should reduce by area ratio
        sub_val = super_val*subFace.volume/superFace.volume
        setfield!(subFace, property, sub_val)
    else
        super_val = getfield(superFace, property)
        sub_val = getfield(subFace, property)
    
        if isa(super_val, Vector)
            # Spectral - copy entire vector
            sub_val .= super_val
        else
            # Grey - direct assignment
            setfield!(subFace, property, super_val)
        end
    end
end