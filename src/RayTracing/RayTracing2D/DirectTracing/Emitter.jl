const WALL_TOLERANCE = 1e-10 # Tolerance for wall detection
const xVecGlobal = SVector(1.0, 0.0)
const yVecGlobal = SVector(0.0, 1.0)

mutable struct Emitter
    type::Symbol  # :surface or :volume
    coarse_index::Int
    fine_index::Int
    wall_index::Int
    energy::G where {G}
end