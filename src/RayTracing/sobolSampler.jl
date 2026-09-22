mutable struct SobolEmitterRNG <: Random.AbstractRNG
    seq::SobolSeq              # the emitter's Sobol sequence (keeps its position across chunks)
    point::Vector{Float64}     # current Sobol point, before the digital shift
    shift::Vector{UInt64}      # digital shift, one 64-bit mask per coordinate
    d::Int                     # numbers consumed per ray
    j::Int                     # coordinates handed out for the current ray
    n::Int                     # rays started so far
end

function SobolEmitterRNG(d::Integer, seed::Integer, bin::Integer, emitter::Integer)
    d >= 1 || error("SobolEmitterRNG: the number of coordinates per ray must be ≥ 1, got $d")
    shift = [_sobolShift(seed, bin, emitter, j) for j in 1:d]   # one mask per coordinate
    # j starts at d ("all coordinates used"), so drawing before nextRay! is an error
    return SobolEmitterRNG(SobolSeq(Int(d)), zeros(Float64, d), shift, Int(d), Int(d), 0)
end

function nextRay!(q::SobolEmitterRNG)
    if q.n == 0
        fill!(q.point, 0.0)             # Sobol index 0; Sobol.jl skips it, the net needs it
    else
        Sobol.next!(q.seq, q.point)     # Sobol index n
    end
    q.n += 1                            # one more ray started
    q.j = 0                             # no coordinates handed out yet
    return nothing
end

nextRay!(::Random.AbstractRNG) = nothing

# Next coordinate of the current point as 64 shifted bits
@inline function _sobolBits(q::SobolEmitterRNG)
    q.j < q.d || error("SobolEmitterRNG: a ray asked for coordinate $(q.j + 1) of a $(q.d)-dimensional " *
                       "point (either nextRay! was not called, or the draw order contract was broken)")
    q.j += 1
    u32 = unsafe_trunc(UInt64, q.point[q.j] * 4294967296.0)   # exact: Sobol points are k / 2^32
    return (u32 << 32) ⊻ q.shift[q.j]                         # low 32 bits are fixed random fill
end

# The tracers draw with rand(rng, G), G being the floating point type of the mesh
Base.rand(q::SobolEmitterRNG, ::Type{Float64}) = (_sobolBits(q) >> 11) * 2.0^-53             # top 53 bits, in [0, 1)
Base.rand(q::SobolEmitterRNG, ::Type{Float32}) = Float32(_sobolBits(q) >> 40) * 2.0f0^-24    # top 24 bits, in [0, 1)
Base.rand(q::SobolEmitterRNG, ::Type{G}) where {G<:AbstractFloat} = G(rand(q, Float64))      # wider types

function sobolDims2D(is_surface::Bool, n_vertices::Integer, samples_depth::Bool)
    d_emit = is_surface ? 3 : (n_vertices == 4 ? 5 : 4)   # wall 3, quad 5, triangle 4
    return d_emit + (samples_depth ? 1 : 0)               # one more for the absorption depth
end

function sobolEmitterRNGs2D(rtm, all_emitters, num_surfaces::Integer, samples_depth::Bool,
                            seed::Integer, bin::Integer)
    ray_rngs = Vector{SobolEmitterRNG}(undef, length(all_emitters))
    for (e, (emitter_key, global_idx)) in enumerate(all_emitters)
        is_surface = global_idx <= num_surfaces                                        # walls come first
        n_vertices = length(rtm.fine_mesh[emitter_key[1]][emitter_key[2]].vertices)    # 3 or 4
        d = sobolDims2D(is_surface, n_vertices, samples_depth)                         # numbers per ray
        ray_rngs[e] = SobolEmitterRNG(d, seed, bin, global_idx)                        # shift unique to this emitter
    end
    return ray_rngs
end

# SplitMix64 finaliser: a stateless, well-mixing 64-bit hash
function _splitmix64(x::UInt64)
    z = x + 0x9e3779b97f4a7c15
    z = (z ⊻ (z >> 30)) * 0xbf58476d1ce4e5b9
    z = (z ⊻ (z >> 27)) * 0x94d049bb133111eb
    return z ⊻ (z >> 31)
end

# Digital shift for one coordinate; chained so that swapping two inputs gives a different mask
function _sobolShift(seed::Integer, bin::Integer, emitter::Integer, coordinate::Integer)
    h = _splitmix64(seed % UInt64)
    h = _splitmix64(h ⊻ (bin % UInt64))
    h = _splitmix64(h ⊻ (emitter % UInt64))
    return _splitmix64(h ⊻ (coordinate % UInt64))
end