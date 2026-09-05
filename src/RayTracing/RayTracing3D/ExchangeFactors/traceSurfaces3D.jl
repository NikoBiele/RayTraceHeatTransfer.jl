# ---------------------------------------------------------------- sampling

@inline function sampleTriangle(a::Point3{G}, b::Point3{G}, c::Point3{G}, rng) where {G}
    s = sqrt(rand(rng, G))
    u = rand(rng, G)
    return (one(G) - s) * a + s * (one(G) - u) * b + s * u * c
end

@inline function sampleFacetPoint(f::Facet3D{G}, rng) where {G}
    if f.nv == 3 || rand(rng, G) < f.tri_split
        return sampleTriangle(f.v[1], f.v[2], f.v[3], rng)
    else
        return sampleTriangle(f.v[1], f.v[3], f.v[4], rng)
    end
end

@inline function sampleCosineDirection(f::Facet3D{G}, rng) where {G}
    r1    = rand(rng, G)
    r2    = rand(rng, G)
    sinth = sqrt(r1)
    costh = sqrt(one(G) - r1)
    phi   = 2 * G(pi) * r2
    p1 = f.e1 * (sinth * cos(phi)) +
           f.e2 * (sinth * sin(phi)) +
           f.inwardNormal * costh
    return Point3{G}(p1[1], p1[2], p1[3])
end

@inline function emitSurfaceRay3D(f::Facet3D{G}, nudge::G, rng) where {G}
    p_emit   = sampleFacetPoint(f, rng) + f.inwardNormal * nudge
    dir_emit = sampleCosineDirection(f, rng)
    return p_emit, dir_emit
end

# ---------------------------------------------------------------- tracing

"""
    traceSurfaces3D(domain, rays_tot; seed, verbose)

First-hit Monte Carlo trace over a transparent 3D surface enclosure,
producing a sparse wavelength-independent `F_raw`. No reflections are
traced: emissivity enters downstream in the GERT solve. Rays that escape
the enclosure are reported by `row_normalize!` as the row-sum deficit.
"""
function traceSurfaces3D(domain::RayTracingDomain3D_surfaces{G,P}, rays_tot::Integer,
                        trace_nudge::G; 
                        nthreads::S=Threads.nthreads(),
                        seeds::Union{Vector{K},K,UnitRange{K},Nothing}=nothing,
                        rngs::Union{Vector{<:Random.AbstractRNG},<:Random.AbstractRNG,Nothing}=nothing,
                        verbose::Bool = true) where {G,S<:Integer,K<:Integer,P<:Integer}

    if nthreads > Threads.nthreads()
        @warn "The number of input threads is higher than the available number of the session."
    end

    if seeds === nothing
        seeds = 1:nthreads
    elseif seeds isa Integer
        seeds = seeds:(seeds + nthreads - 1)
    end
    length(seeds) == nthreads ||
        error("got $(length(seeds)) seeds for $nthreads threads; supply one per thread, a single starting seed, or none")
    length(unique(seeds)) == nthreads ||
        error("seeds must be distinct; threads sharing a seed trace identical rays")

    if rngs === nothing
        rngs = [Random.Xoshiro(s) for s in seeds]
    elseif rngs isa Random.AbstractRNG
        rngs = [deepcopy(rngs) for _ in 1:nthreads]
    end
    length(rngs) == nthreads ||
        error("got $(length(rngs)) rngs for $nthreads threads")
    length(unique(objectid.(rngs))) == nthreads ||
        error("rngs must be distinct objects; `fill(Xoshiro(1), n)` aliases one RNG across all threads")

    domain.flattened || flatten!(domain)

    facets = domain.facets
    tris   = domain.tris
    bvh    = domain.bvh
    nudge  = trace_nudge

    num_emitters     = length(facets)
    rays_per_emitter = max(1, cld(Int(rays_tot), num_emitters))
    inv_rays         = one(G) / G(rays_per_emitter)

    verbose && println("  Using $nthreads threads, " *
                       "$rays_per_emitter rays from each of $num_emitters elements")

    emitters_per_thread = div(num_emitters, nthreads)
    remainder           = num_emitters % nthreads
    thread_assignments  = Vector{UnitRange{Int}}(undef, nthreads)
    start_idx = 1
    for tid in 1:nthreads
        thread_size = emitters_per_thread + (tid <= remainder ? 1 : 0)
        thread_assignments[tid] = start_idx:(start_idx + thread_size - 1)
        start_idx += thread_size
    end

    Is = [Int[] for _ in 1:nthreads]
    Js = [Int[] for _ in 1:nthreads]
    Vs = [G[]   for _ in 1:nthreads]

    progress  = Progress(num_emitters; dt = 1, desc = "  Tracing progress: ")
    completed = Threads.Atomic{Int}(0)

    @threads for tid in 1:nthreads
        Il, Jl, Vl = Is[tid], Js[tid], Vs[tid]
        row   = Dict{Int,Int}()          # absorber → count, reused per emitter
        stack = Int32[]                  # BVH traversal stack, reused per emitter
        sizehint!(stack, 64)
        local_rng   = rngs[tid]
        Random.seed!(rngs[tid], seeds[tid])

        for global_idx in thread_assignments[tid]
            f = facets[global_idx]
            empty!(row)

            for _ in 1:rays_per_emitter
                p_emit, dir_emit = emitSurfaceRay3D(f, nudge, local_rng)
                a, _ = closestHit(bvh, tris, p_emit, dir_emit, stack)
                a < 0 && continue                       # escaped the enclosure
                row[a] = get(row, a, 0) + 1
            end

            for (j, c) in row
                push!(Il, Int(global_idx)); push!(Jl, Int(j)); push!(Vl, c * inv_rays)
            end

            Threads.atomic_add!(completed, 1)
            tid == 1 && update!(progress, completed[])
        end
    end
    finish!(progress)

    F_raw = sparse(reduce(vcat, Is), reduce(vcat, Js), reduce(vcat, Vs),
                   num_emitters, num_emitters)
    empty!.(Is); empty!.(Js); empty!.(Vs)
    GC.gc()

    domain.F_raw = row_normalize!(F_raw, rays_per_emitter)
    return domain
end

@inline function intersectTri(t::Tri3D{G}, o::Point3{G}, d::Point3{G}, tmax::G) where {G}
    pv  = cross(d, t.ac)
    det = dot(t.ab, pv)
    abs(det) < 1e-14 && return (false, tmax)          # ray parallel to triangle
    inv = one(G) / det
    tv  = o - t.a
    u   = dot(tv, pv) * inv
    (u < 0 || u > 1) && return (false, tmax)
    qv = cross(tv, t.ab)
    v  = dot(d, qv) * inv
    (v < 0 || u + v > 1) && return (false, tmax)
    tt = dot(t.ac, qv) * inv
    (tt <= 0 || tt >= tmax) && return (false, tmax)
    return (true, tt)
end

@inline function slabHit(node::BVHNode{G}, o::Point3{G}, invd::Point3{G}, tmax::G) where {G}
    tmin = zero(G); tmx = tmax
    @inbounds for ax in 1:3
        t0 = (node.bmin[ax] - o[ax]) * invd[ax]
        t1 = (node.bmax[ax] - o[ax]) * invd[ax]
        t0 > t1 && ((t0, t1) = (t1, t0))
        tmin = max(tmin, t0); tmx = min(tmx, t1)
        tmin > tmx && return false
    end
    return true
end

function closestHit(bvh::BVH{G}, tris::Vector{Tri3D{G}},
                    o::Point3{G}, d::Point3{G},
                    stack::Vector{Int32}; tmax::G = G(Inf)) where {G}

    invd  = Point3{G}(one(G)/d[1], one(G)/d[2], one(G)/d[3])
    bestT = tmax; bestElem = Int32(-1)

    empty!(stack); push!(stack, Int32(1))
    @inbounds while !isempty(stack)
        node = bvh.nodes[pop!(stack)]
        slabHit(node, o, invd, bestT) || continue
        if node.left < 0
            for k in node.first : node.first + node.count - 1
                ti = bvh.order[k]
                hit, tt = intersectTri(tris[ti], o, d, bestT)
                if hit
                    bestT = tt; bestElem = tris[ti].element
                end
            end
        else
            push!(stack, node.left)
            push!(stack, node.left + Int32(1))
        end
    end
    return bestElem, bestT
end