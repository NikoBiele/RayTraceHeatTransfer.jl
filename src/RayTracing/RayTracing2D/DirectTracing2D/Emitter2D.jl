struct Emitter2D{G}
    type::Symbol  # :surface or :volume
    coarse_index::Int
    fine_index::Int
    wall_index::Int
    energy::G     # emitted power of this element in the traced bin [W]
end