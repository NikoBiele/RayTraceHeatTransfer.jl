# Grey directional solver.
#
# State J[k, b]: power leaving element k in direction bin b (N × A, stored bin-major).
# With G the angular exchange factors, S[k, b] = Σ_j G[b][k, j] the emission shares
# and F_b = G[b] ./ S[:, b] the row-stochastic transport of bin b:
#
#     incident   g[k, b]  = Σ_j F_b[j, k] J[j, b]
#     J[k, b′] − Σ_b C_k[b, b′] g[k, b] = h0[k] · S[k, b′]
#
#     known T:  C_k = ρ_k Φ_k                       h0 = emission
#     known q:  C_k = (1 − ρ_k)·𝟙S_kᵀ + ρ_k Φ_k     h0 = source
#
# Φ_k = d_k·𝟙S_kᵀ + t_k·Table_k: a diffuse part (rank 1, applied in O(A)) plus an
# optional table shared between all elements with the same descriptor. With every
# Φ_k diffuse, summing over b′ gives exactly the grey system I − diag(coeff)·Fᵀ.

struct DirectionalSystemOp
    N::Int
    A::Int
    G::Vector{SparseMatrixCSC{Float64,Int}}
    S::Matrix{Float64}            # N × A emission shares (row sums of G)
    invS::Matrix{Float64}         # 1 ./ S, 0 where S == 0
    refl::Vector{Float64}         # ρ_k: wall reflectivity or scattering albedo
    reemit::Vector{Float64}       # 1 − ρ_k for flux-specified elements, else 0
    diff_w::Vector{Float64}       # d_k
    tab_w::Vector{Float64}        # t_k
    tab_idx::Vector{Int}          # index into tables, 0 = none
    tables::Vector{Matrix{Float64}}
    g::Matrix{Float64}            # scratch: incident power N × A
    xs::Matrix{Float64}           # scratch
    tmp::Vector{Float64}          # scratch for the 5-argument mul!
end

Base.size(op::DirectionalSystemOp)            = (op.N * op.A, op.N * op.A)
Base.size(op::DirectionalSystemOp, ::Integer) = op.N * op.A
Base.eltype(::DirectionalSystemOp)            = Float64

# op.g ← incident power for the state X (N × A)
function _incident!(op::DirectionalSystemOp, X::AbstractMatrix)
    @. op.xs = X * op.invS
    Threads.@threads for b in 1:op.A
        mul!(view(op.g, :, b), transpose(op.G[b]), view(op.xs, :, b))
    end
    return op.g
end

# out ← redistribution of the incident power held in op.g
function _redistribute!(out::AbstractMatrix, op::DirectionalSystemOp)
    g, S, A = op.g, op.S, op.A
    Threads.@threads for k in 1:op.N
        gt = 0.0
        @inbounds for b in 1:A
            gt += g[k, b]
        end
        cd = (op.refl[k] * op.diff_w[k] + op.reemit[k]) * gt
        @inbounds for bp in 1:A
            out[k, bp] = cd * S[k, bp]
        end
        ti = op.tab_idx[k]
        ct = op.refl[k] * op.tab_w[k]
        if ti != 0 && ct != 0.0
            Φ = op.tables[ti]
            @inbounds for bp in 1:A
                s = 0.0
                for b in 1:A
                    s += Φ[b, bp] * g[k, b]
                end
                out[k, bp] += ct * s
            end
        end
    end
    return out
end

function LinearAlgebra.mul!(y::AbstractVector, op::DirectionalSystemOp, x::AbstractVector)
    X = reshape(x, op.N, op.A)
    Y = reshape(y, op.N, op.A)
    _incident!(op, X)
    _redistribute!(Y, op)
    @. y = x - y
    return y
end

function LinearAlgebra.mul!(y::AbstractVector, op::DirectionalSystemOp, x::AbstractVector,
                            α::Number, β::Number)
    mul!(op.tmp, op, x)
    @. y = α * op.tmp + β * y
    return y
end

# descriptor → (diffuse weight, table weight, cache key, table builder)
_element_redistribution(dm, ::IsotropicScattering, φn) = (1.0, 0.0, nothing, nothing)
_element_redistribution(dm, ::DiffuseReflection, φn)   = (1.0, 0.0, nothing, nothing)
_element_redistribution(dm, d::HenyeyGreenstein, φn) =
    (0.0, 1.0, (:hg, d.g), () -> redistribution_table(dm, d))
_element_redistribution(dm, d::TabulatedScattering, φn) =
    (0.0, 1.0, (:tabulated_scattering, objectid(d.table)), () -> redistribution_table(dm, d))
function _element_redistribution(dm, d::SpecularReflection, φn)
    d.specularity == 0 && return (1.0, 0.0, nothing, nothing)
    return (1.0 - d.specularity, d.specularity, (:mirror, round(φn; digits = 12)), () -> mirror_table(dm, φn))
end
_element_redistribution(dm, d::TabulatedReflection, φn) =
    (0.0, 1.0, (:tabulated_reflection, objectid(d.table), round(φn; digits = 12)),
     () -> redistribution_table(dm, d, φn))

# Operator, right-hand-side scale h0, workspace and albedo for one (grey or single spectral bin) system
function _directional_setup(rtm::RayTracingDomain2D, G::Vector{SparseMatrixCSC{Float64,Int}};
                            spectral_bin::Int = 1)
    dm = rtm.directional_model
    validate(dm)
    A  = n_angular_bins(dm)
    length(G) == A || error("G has $(length(G)) direction bins but the directional model has $A; retrace")
    ns = length(rtm.surface_mapping)
    nv = length(rtm.volume_mapping)
    N  = ns + nv

    # boundary data exactly as in the grey solver
    ws = OneMatrixWorkspace{Float64}(ns, nv)
    E_known = ws.work_vec1
    Q_known = ws.work_vec2
    populateWorkspace!(ws, rtm, spectral_bin)
    computeEmissivePowersVariable!(rtm, ws, E_known, Q_known)
    b = zeros(N)
    for s in 1:ns
        b[s] = 1.0 - ws.epsw[s]
    end
    for v in 1:nv
        b[ns + v] = ws.omega_g[v]
    end
    known_q = vcat(ws.bin_Qw_known, ws.bin_Qg_known) .== 1
    h0      = [known_q[i] ? Q_known[i] : E_known[i] for i in 1:N]
    reemit  = [known_q[i] ? 1.0 - b[i] : 0.0 for i in 1:N]

    # shares of the G in use
    S = zeros(N, A)
    for a in 1:A
        S[:, a] .= vec(sum(G[a], dims = 2))
    end
    invS = [S[i, a] > 0 ? 1.0 / S[i, a] : 0.0 for i in 1:N, a in 1:A]

    # per-element redistribution; tables shared between equal descriptors
    diff_w = ones(N); tab_w = zeros(N); tab_idx = zeros(Int, N)
    tables = Matrix{Float64}[]
    cache  = Dict{Any,Int}()
    function assign!(k, d, φn)
        d isa Vector && (d = d[spectral_bin])
        dw, tw, key, build = _element_redistribution(dm, d, φn)
        diff_w[k] = dw; tab_w[k] = tw
        if key !== nothing
            idx = get(cache, key, 0)
            if idx == 0
                push!(tables, build())
                idx = length(tables)
                cache[key] = idx
            end
            tab_idx[k] = idx
        end
    end
    for ((ci, fi, wi), s) in rtm.surface_mapping
        f = rtm.fine_mesh[ci][fi]
        assign!(s, f.reflection[wi], wall_normal_azimuth(f, wi))
    end
    for ((ci, fi), v) in rtm.volume_mapping
        assign!(ns + v, rtm.fine_mesh[ci][fi].phase, 0.0)
    end

    op = DirectionalSystemOp(N, A, G, S, invS, b, reemit, diff_w, tab_w, tab_idx, tables,
                             zeros(N, A), zeros(N, A), zeros(N * A))
    return op, h0, ws, b, known_q
end

"""
    equilibriumGreyDirectional2D!(rtm, F; verbose = true)

Grey solve with directional redistribution (`rtm.directional_model` set): the
faces' `phase` and `reflection` descriptors decide how scattered and reflected
power is distributed over the direction bins. Uses `rtm.G_raw` when `F` is
`rtm.F_raw` itself and `rtm.G_smooth` otherwise. Element totals are written to
the faces as by the grey solver; the angular solution is stored in
`rtm.J` (N × A, power leaving each element per direction bin).
"""
function equilibriumGreyDirectional2D!(rtm::RayTracingDomain2D, F::AbstractMatrix; verbose::Bool = true)
    verbose && println("=== Grey directional solver ===")
    G = F === rtm.F_raw ? rtm.G_raw : rtm.G_smooth
    G === nothing && error(F === rtm.F_raw ?
        "no angular exchange factors: trace with domain(N_rays; method = :pathlength) after setting directional_model" :
        "G_smooth not computed: call smooth!(domain), or pass domain.F_raw itself to solve on the raw exchange factors")
    G isa Vector{SparseMatrixCSC{Float64,Int}} ||
        error("spectral angular exchange factors on a grey solve; retrace the domain")

    op, h0, ws, b = _directional_setup(rtm, G)
    N, A = op.N, op.A
    verbose && println("  $N elements × $A direction bins = $(N * A) unknowns, ",
                       "$(length(op.tables)) shared redistribution table(s)")

    h   = vec(h0 .* op.S)
    wkr = GmresWorkspace(op, h; memory = 50)
    gmres!(wkr, op, h; restart = true, rtol = 1e-12, atol = 0.0)
    wkr.stats.solved || @warn "directional GMRES did not reach rtol = 1e-12" wkr.stats.niter
    verbose && println("  GMRES iterations: $(wkr.stats.niter)")

    J = reshape(copy(wkr.x), N, A)
    _incident!(op, J)
    j    = vec(sum(J, dims = 2))
    gtot = vec(sum(op.g, dims = 2))
    r    = b .* gtot
    Abs  = (1.0 .- b) .* gtot

    T = zeros(N)
    computeTemperaturesVariable!(rtm, ws, T, j, r)
    writeResultsToDomain!(rtm, j, Abs, r; T = T, spectral_bin = 1)
    rtm.J = J
    rtm.energy_error  = sum(j - r - Abs) / (sum(j) > 1000 * eps(Float64) ? sum(j) : one(Float64))

    if verbose
        show(stdout, MIME"text/plain"(), rtm)
        println()
    end
    return nothing
end