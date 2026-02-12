mutable struct Emitter
    type::Symbol  # :surface or :volume
    coarse_index::Int
    fine_index::Int
    wall_index::Int
    energy::G where {G}
end