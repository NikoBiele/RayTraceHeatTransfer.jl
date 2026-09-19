# Directional (angular) models.
#
# A directional domain carries one `directional_model`, which fixes the angular
# bins. Element-level descriptors (AbstractPhaseFunction for volumes,
# AbstractWallReflection for walls) say how incident power is redistributed over
# those bins; the tables are built here at the model's resolution.
#
# Direction convention (2D tracer): a ray direction is the in-plane projection
# (d₁, d₂) of a 3D unit vector; azimuth φ = atan(d₂, d₁), out-of-plane
# |μ_z| = √(1 − d₁² − d₂²). Bin index a = (iz − 1)·n_azimuth + iφ.
#
# A redistribution table Φ[b, b′] (incident bin b → outgoing bin b′) is valid when
#     Φ ≥ 0,   Σ_b′ Φ[b, b′] = 1                       (conservation)
#     Σ_b p[rev(b)] Φ[b, b′] = p[b′]                   (detailed balance)
# with p the element's emission shares. With reciprocal exchange factors the
# second condition makes an isothermal enclosure exact for any table.

n_angular_bins(m::AngularBins) = m.n_azimuth * m.n_polar

function validate(m::AngularBins)
    (m.n_azimuth >= 2 && iseven(m.n_azimuth)) ||
        error("AngularBins: n_azimuth must be even and ≥ 2 (every direction bin needs its reverse), got $(m.n_azimuth)")
    m.n_polar >= 1 || error("AngularBins: n_polar must be ≥ 1, got $(m.n_polar)")
    return nothing
end

describe(m::AngularBins) =
    "AngularBins: $(m.n_azimuth) azimuthal × $(m.n_polar) out-of-plane = $(n_angular_bins(m)) direction bins"

# bin of a stored ray direction (any indexable (d₁, d₂))
function angular_bin(m::AngularBins, d)
    d1, d2 = Float64(d[1]), Float64(d[2])
    φ  = atan(d2, d1)
    μz = sqrt(max(0.0, 1.0 - d1^2 - d2^2))
    iφ = clamp(floor(Int, (φ + π) / (2π) * m.n_azimuth) + 1, 1, m.n_azimuth)
    iz = clamp(floor(Int, μz * m.n_polar) + 1, 1, m.n_polar)
    return (iz - 1) * m.n_azimuth + iφ
end

# bin of the reversed direction: φ → φ + π, |μ_z| unchanged
function reversed_bin(m::AngularBins, a::Int)
    iz = (a - 1) ÷ m.n_azimuth + 1
    iφ = (a - 1) % m.n_azimuth + 1
    return (iz - 1) * m.n_azimuth + mod1(iφ + m.n_azimuth ÷ 2, m.n_azimuth)
end

azimuth_edges(m::AngularBins, iφ::Int) =
    (-π + (iφ - 1) * 2π / m.n_azimuth, -π + iφ * 2π / m.n_azimuth)

# ∫ max(cos(φ − φ0), 0) dφ over [φlo, φhi]  (|φhi − φlo| ≤ 2π)
function _cos_lobe_integral(φlo::Real, φhi::Real, φ0::Real)
    c = mod(φ0 + π, 2π) - π                              # lobe centre in [−π, π)
    total = 0.0
    for k in -2:2
        lo = max(φlo - c + 2π * k, -π / 2)
        hi = min(φhi - c + 2π * k,  π / 2)
        hi > lo && (total += sin(hi) - sin(lo))
    end
    return total
end

# ∫ √(1 − μ²) dμ
_sqrt_one_minus_sq_primitive(μ::Real) = 0.5 * (μ * sqrt(max(0.0, 1.0 - μ^2)) + asin(μ))

# ---------------------------------------------------------- emission shares ---
# Fraction of an element's emission in each bin. Volumes: isotropic.
emission_shares(m::AngularBins) = fill(1.0 / n_angular_bins(m), n_angular_bins(m))

# Walls: Lambertian about the in-plane emission normal at azimuth φn,
# pdf = (d·n)/π with d·n = √(1 − μ_z²)·cos(φ − φn). Bins behind the wall get exactly 0.
function emission_shares(m::AngularBins, φn::Real)
    p = zeros(n_angular_bins(m))
    for iz in 1:m.n_polar, iφ in 1:m.n_azimuth
        φlo, φhi = azimuth_edges(m, iφ)
        Iμ = _sqrt_one_minus_sq_primitive(iz / m.n_polar) - _sqrt_one_minus_sq_primitive((iz - 1) / m.n_polar)
        p[(iz - 1) * m.n_azimuth + iφ] = (2 / π) * Iμ * _cos_lobe_integral(φlo, φhi, φn)
    end
    p[p .< 1e-12] .= 0.0
    p ./= sum(p)
    return p
end

# ------------------------------------------------------------ table checks ---
function check_redistribution_table(Φ::AbstractMatrix, p::AbstractVector, m::AngularBins; atol::Real = 1e-10)
    A = n_angular_bins(m)
    size(Φ) == (A, A) || error("redistribution table is $(size(Φ)) but the directional model has $A bins")
    length(p) == A || error("emission shares have $(length(p)) entries but the directional model has $A bins")
    all(>=(0), Φ) || error("redistribution table has negative entries")
    rowerr = maximum(abs.(vec(sum(Φ, dims = 2)) .- 1.0))
    rowerr <= atol || error("redistribution table rows must sum to 1 (max deviation $rowerr)")
    for bp in 1:A
        s = 0.0
        for b in 1:A
            s += p[reversed_bin(m, b)] * Φ[b, bp]
        end
        abs(s - p[bp]) <= atol ||
            error("redistribution table violates detailed balance in bin $bp: Σ_b p[rev b]·Φ[b, b′] = $s, p[b′] = $(p[bp])")
    end
    return nothing
end

# ------------------------------------------------------------------ volumes ---
redistribution_table(m::AngularBins, ::IsotropicScattering) =
    fill(1.0 / n_angular_bins(m), n_angular_bins(m), n_angular_bins(m))

# Bin-averaged Henyey–Greenstein kernel (nq × nq sub-samples per bin, ±z folded),
# then symmetric scaling to a symmetric doubly stochastic table: rows give
# conservation, columns give detailed balance against the isotropic shares.
function redistribution_table(m::AngularBins, pf::HenyeyGreenstein; nq::Int = 4)
    Aφ, Az = m.n_azimuth, m.n_polar
    A, nQ, g = Aφ * Az, nq * nq, pf.g
    dirs = Matrix{NTuple{3,Float64}}(undef, A, nQ)
    for iz in 1:Az, iφ in 1:Aφ
        a = (iz - 1) * Aφ + iφ
        q = 0
        for qz in 1:nq, qφ in 1:nq
            q += 1
            φ = -π + (iφ - 1 + (qφ - 0.5) / nq) * 2π / Aφ
            μ = (iz - 1 + (qz - 0.5) / nq) / Az
            s = sqrt(1.0 - μ^2)
            dirs[a, q] = (s * cos(φ), s * sin(φ), μ)
        end
    end
    Φ = zeros(A, A)
    for a in 1:A, b in 1:A
        acc = 0.0
        for qa in 1:nQ, qb in 1:nQ
            da, db = dirs[a, qa], dirs[b, qb]
            inplane = da[1] * db[1] + da[2] * db[2]
            for sgn in (1.0, -1.0)
                c = clamp(inplane + sgn * da[3] * db[3], -1.0, 1.0)
                acc += (1.0 - g^2) / (1.0 + g^2 - 2.0 * g * c)^1.5
            end
        end
        Φ[a, b] = acc
    end
    for _ in 1:1000
        r = vec(sum(Φ, dims = 2))
        maximum(abs.(r .- 1.0)) < 1e-14 && break
        d = 1.0 ./ sqrt.(r)
        Φ .= d .* Φ .* d'
    end
    check_redistribution_table(Φ, emission_shares(m), m)
    return Φ
end

function redistribution_table(m::AngularBins, pf::TabulatedScattering)
    check_redistribution_table(pf.table, emission_shares(m), m)
    return copy(pf.table)
end

# -------------------------------------------------------------------- walls ---
function redistribution_table(m::AngularBins, ::DiffuseReflection, φn::Real)
    p = emission_shares(m, φn)
    return repeat(reshape(p, 1, :), length(p), 1)          # Φ[b, b′] = p[b′]
end

# Ideal mirror about a wall with emission-normal azimuth φn: φ ↦ 2φn + π − φ,
# |μ_z| kept. Entry [b, b′] is the share of the incident flux in bin b (weight
# (−d·n)⁺, i.e. a cosine lobe about φn + π) whose mirror image falls in b′ —
# computed exactly from interval overlaps, for any wall orientation. It reduces
# to a permutation when the bin edges lie on the wall tangent, and satisfies
# detailed balance against the Lambertian shares by construction. Bins carrying
# no incident flux get the diffuse row (never used, keeps every row stochastic).
function mirror_table(m::AngularBins, φn::Real)
    Aφ, Az = m.n_azimuth, m.n_polar
    φin = φn + π
    Pφ = zeros(Aφ, Aφ)
    incident = falses(Aφ)
    for iφ in 1:Aφ
        lo, hi = azimuth_edges(m, iφ)
        q = _cos_lobe_integral(lo, hi, φin)
        q > 1e-12 || continue
        incident[iφ] = true
        for iφp in 1:Aφ
            lop, hip = azimuth_edges(m, iφp)
            plo, phi = 2φn + π - hip, 2φn + π - lop        # preimage of bin iφp under the mirror
            acc = 0.0
            for k in -2:2
                a = max(lo, plo + 2π * k)
                b = min(hi, phi + 2π * k)
                b > a && (acc += _cos_lobe_integral(a, b, φin))
            end
            Pφ[iφ, iφp] = acc / q
        end
        row = view(Pφ, iφ, :)
        row[row .< 1e-12] .= 0.0
        row ./= sum(row)
    end
    p = emission_shares(m, φn)
    P = zeros(Aφ * Az, Aφ * Az)
    for iz in 1:Az, iφ in 1:Aφ
        b = (iz - 1) * Aφ + iφ
        if incident[iφ]
            for iφp in 1:Aφ
                P[b, (iz - 1) * Aφ + iφp] = Pφ[iφ, iφp]
            end
        else
            P[b, :] .= p
        end
    end
    return P
end

function redistribution_table(m::AngularBins, r::SpecularReflection, φn::Real)
    p = emission_shares(m, φn)
    Φ = r.specularity .* mirror_table(m, φn) .+ (1.0 - r.specularity) .* repeat(reshape(p, 1, :), length(p), 1)
    check_redistribution_table(Φ, p, m)
    return Φ
end

function redistribution_table(m::AngularBins, r::TabulatedReflection, φn::Real)
    check_redistribution_table(r.table, emission_shares(m, φn), m)
    return copy(r.table)
end


# ------------------------------------------------------ domain-level shares ---
# azimuth of the emission normal of wall `wi` of a fine face: the normal used by
# emitSurfaceRay2D is (−t₂, t₁) for the wall tangent t = p₂ − p₁
function wall_normal_azimuth(f::PolyVolume2D, wi::Int)
    p1 = f.vertices[wi]
    p2 = f.vertices[mod1(wi + 1, length(f.vertices))]
    t1, t2 = p2[1] - p1[1], p2[2] - p1[2]
    return atan(t1, -t2)
end

# N × A matrix of analytic emission shares in global element order
# (walls Lambertian about their normal, volumes isotropic)
function emission_share_matrix(rtm::RayTracingDomain2D)
    dm = rtm.directional_model
    dm === nothing && error("the domain has no directional_model")
    ns = length(rtm.surface_mapping)
    N  = ns + length(rtm.volume_mapping)
    p  = zeros(N, n_angular_bins(dm))
    for ((ci, fi, wi), s) in rtm.surface_mapping
        p[s, :] .= emission_shares(dm, wall_normal_azimuth(rtm.fine_mesh[ci][fi], wi))
    end
    pv = emission_shares(dm)
    for (_, v) in rtm.volume_mapping
        p[ns + v, :] .= pv
    end
    return p
end