@inline bmin2(a::Point3{G}, b::Point3{G}) where {G} =
    Point3{G}(min(a[1],b[1]), min(a[2],b[2]), min(a[3],b[3]))
@inline bmax2(a::Point3{G}, b::Point3{G}) where {G} =
    Point3{G}(max(a[1],b[1]), max(a[2],b[2]), max(a[3],b[3]))

@inline function boxArea(lo::Point3{G}, hi::Point3{G}) where {G}
    d  = hi - lo
    dx = max(d[1], zero(G)); dy = max(d[2], zero(G)); dz = max(d[3], zero(G))
    return 2*(dx*dy + dy*dz + dz*dx)
end

@inline function triBounds(t::Tri3D{G}) where {G}
    p1 = t.a; p2 = t.a + t.ab; p3 = t.a + t.ac
    lo = bmin2(bmin2(p1, p2), p3)
    hi = bmax2(bmax2(p1, p2), p3)
    return lo, hi
end

emptyNode(::Type{G}) where {G} =
    BVHNode{G}(Point3{G}(0,0,0), Point3{G}(0,0,0), Int32(0), Int32(0), Int32(0))

const BVH_BINS     = 12
const BVH_MAX_LEAF = 8

function buildBVH(tris::Vector{Tri3D{G}}) where {G}
    n = length(tris)
    n == 0 && error("buildBVH: no primitives — call flatten! first")

    lo  = Vector{Point3{G}}(undef, n)
    hi  = Vector{Point3{G}}(undef, n)
    cen = Vector{Point3{G}}(undef, n)
    for i in 1:n
        lo[i], hi[i] = triBounds(tris[i])
        cen[i] = (lo[i] + hi[i]) / 2
    end

    order = Int32.(collect(1:n))
    nodes = BVHNode{G}[emptyNode(G)]          # root placeholder
    buildNode!(nodes, order, lo, hi, cen, 1, 1, n)
    return BVH{G}(nodes, order)
end

function buildNode!(nodes::Vector{BVHNode{G}}, order, lo, hi, cen,
                    nodeIdx::Int, first::Int, count::Int) where {G}

    blo = Point3{G}(Inf,Inf,Inf);    bhi = Point3{G}(-Inf,-Inf,-Inf)
    clo = Point3{G}(Inf,Inf,Inf);    chi = Point3{G}(-Inf,-Inf,-Inf)
    for k in first:first+count-1
        i = order[k]
        blo = bmin2(blo, lo[i]);  bhi = bmax2(bhi, hi[i])
        clo = bmin2(clo, cen[i]); chi = bmax2(chi, cen[i])
    end

    makeLeaf() = (nodes[nodeIdx] =
        BVHNode{G}(blo, bhi, Int32(-1), Int32(first), Int32(count)))

    count <= 2 && return makeLeaf()

    ext  = chi - clo
    axis = argmax((ext[1], ext[2], ext[3]))
    ext[axis] <= 0 && return makeLeaf()       # coincident centroids

    scale = BVH_BINS / ext[axis]
    binof(i) = clamp(Int(floor((cen[i][axis] - clo[axis]) * scale)) + 1, 1, BVH_BINS)

    binLo = fill(Point3{G}(Inf,Inf,Inf),   BVH_BINS)
    binHi = fill(Point3{G}(-Inf,-Inf,-Inf), BVH_BINS)
    binN  = zeros(Int, BVH_BINS)
    for k in first:first+count-1
        i = order[k]; b = binof(i)
        binLo[b] = bmin2(binLo[b], lo[i]); binHi[b] = bmax2(binHi[b], hi[i])
        binN[b] += 1
    end

    leftArea = zeros(G, BVH_BINS); leftN = zeros(Int, BVH_BINS)
    l_lo = Point3{G}(Inf,Inf,Inf); l_hi = Point3{G}(-Inf,-Inf,-Inf); acc = 0
    for b in 1:BVH_BINS
        l_lo = bmin2(l_lo, binLo[b]); l_hi = bmax2(l_hi, binHi[b]); acc += binN[b]
        leftArea[b] = boxArea(l_lo, l_hi); leftN[b] = acc
    end

    bestCost = G(Inf); bestBin = 0
    r_lo = Point3{G}(Inf,Inf,Inf); r_hi = Point3{G}(-Inf,-Inf,-Inf); accR = 0
    for b in BVH_BINS:-1:2
        r_lo = bmin2(r_lo, binLo[b]); r_hi = bmax2(r_hi, binHi[b]); accR += binN[b]
        nL = leftN[b-1]
        (nL == 0 || accR == 0) && continue
        cost = leftArea[b-1]*nL + boxArea(r_lo, r_hi)*accR
        if cost < bestCost
            bestCost = cost; bestBin = b
        end
    end

    if bestBin == 0 || (count <= BVH_MAX_LEAF && bestCost >= boxArea(blo,bhi)*count)
        return makeLeaf()
    end

    mid = first
    for k in first:first+count-1
        if binof(order[k]) < bestBin
            order[k], order[mid] = order[mid], order[k]
            mid += 1
        end
    end
    nLeft = mid - first
    if nLeft == 0 || nLeft == count          # degenerate split → median fallback
        nLeft = count ÷ 2
        mid   = first + nLeft
    end

    leftIdx = length(nodes) + 1
    push!(nodes, emptyNode(G)); push!(nodes, emptyNode(G))
    nodes[nodeIdx] = BVHNode{G}(blo, bhi, Int32(leftIdx), Int32(0), Int32(0))

    buildNode!(nodes, order, lo, hi, cen, leftIdx,     first, nLeft)
    buildNode!(nodes, order, lo, hi, cen, leftIdx + 1, mid,   count - nLeft)
    return nothing
end