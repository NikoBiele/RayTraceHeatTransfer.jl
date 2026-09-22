
# Spectral directional solver — dense path.
#
# Per spectral bin k the state is J_k ∈ ℝ^{N·A} (power leaving each element per
# direction bin, bin-major). The bordered system of the spectral Woodbury solver
# keeps its form with
#
#     𝐃_k  (N·A × N·A)   I − L_k T_k, redistribution with ρ only (no re-emission)
#     𝐌_k  (N × N·A)     border: total outgoing − coeff · total incident
#     P_k  (N·A × N)     spread an element's emission over its direction shares
#     R    (N × N·A)     sum over direction bins
#
# so that block rows read 𝐃_k J_k = P_k e_k and the border Σ_k 𝐌_k J_k = h.
# 𝐃_k does not depend on temperature: it is assembled and factorised once.
# Identities used as tests: R 𝐃_k P_k = I − diag(ρ) F_kᵀ and 𝐌_k P_k = I − diag(coeff) F_kᵀ
# for any valid tables, with F_k = Σₐ G_k[a].

const DIRECTIONAL_DENSE_LIMIT = 6000      # largest N·A handled by the dense path

# Dense 𝐃 and 𝐌 for one spectral bin from the matrix-free operator's data.
# G[b] is CSC: column k holds the senders j of receiver k.
function _directional_dense_blocks(op::DirectionalSystemOp, coeff::AbstractVector)
    N, A = op.N, op.A
    n = N * A
    n <= DIRECTIONAL_DENSE_LIMIT ||
        error("the dense spectral-directional path handles N·A ≤ $DIRECTIONAL_DENSE_LIMIT per spectral bin, " *
              "got $N elements × $A direction bins = $n; use fewer direction bins or a coarser mesh")
    length(coeff) == N || error("coeff has $(length(coeff)) entries for $N elements")

    D = Matrix{Float64}(I, n, n)
    M = zeros(N, n)
    for b in 1:A
        for i in 1:N
            M[i, (b - 1) * N + i] = 1.0
        end
        Gb   = op.G[b]
        rows = rowvals(Gb)
        vals = nonzeros(Gb)
        for k in 1:N
            ρk = op.refl[k]
            dk = op.diff_w[k]
            tk = op.tab_w[k]
            ti = op.tab_idx[k]
            for idx in nzrange(Gb, k)
                j = rows[idx]
                f = vals[idx] * op.invS[j, b]              # F_b[j, k]
                f == 0.0 && continue
                col = (b - 1) * N + j
                M[k, col] -= coeff[k] * f
                ρk == 0.0 && continue
                for bp in 1:A
                    c = dk * op.S[k, bp]
                    (ti != 0 && tk != 0.0) && (c += tk * op.tables[ti][b, bp])
                    D[(bp - 1) * N + k, col] -= ρk * c * f
                end
            end
        end
    end
    return D, M
end

# P e: emission of every element spread over its direction shares
_spread_emission(op::DirectionalSystemOp, e::AbstractVector) = vec(op.S .* e)

# R x: sum over direction bins
_sum_bins(op::DirectionalSystemOp, x::AbstractVector) = vec(sum(reshape(x, op.N, op.A), dims = 2))


"""
    equilibriumSpectralDirectional2D!(rtm, F_matrices; max_iters = 1000, convergence_tol = 1e-12, verbose = true)

Spectral solve with directional redistribution (dense path, N·A ≤ $DIRECTIONAL_DENSE_LIMIT per
spectral bin). The bordered system of the spectral Woodbury solver is kept, with each
bin's block replaced by the directional operator; that operator does not depend on
temperature, so it is LU-factorised once and every sweep costs one solve per bin.
Plain fixed-point sweeps on the emissive power (no Newton acceleration yet).
Uses `rtm.G_raw` when `F_matrices` is `rtm.F_raw` itself and `rtm.G_smooth` otherwise.
Per-bin element totals are written to the faces as by the spectral solver; the angular
solution is stored in `rtm.J` (one N × A matrix per spectral bin).
"""
function equilibriumSpectralDirectional2D!(rtm::RayTracingDomain2D, F_matrices::AbstractVector;
                                           max_iters::Integer = 1000, convergence_tol = 1e-12,
                                           verbose::Bool = true)
    verbose && println("=== Spectral directional solver (dense path) ===")
    validateSpectralSetup(rtm)
    G_all = F_matrices === rtm.F_raw ? rtm.G_raw : rtm.G_smooth
    G_all === nothing && error(F_matrices === rtm.F_raw ?
        "no angular exchange factors: trace with domain(N_rays; method = :pathlength) after setting directional_model" :
        "G_smooth not computed: call smooth!(domain), or pass domain.F_raw itself to solve on the raw exchange factors")
    G_all isa Vector{Vector{SparseMatrixCSC{Float64,Int}}} ||
        error("grey angular exchange factors on a spectral solve; retrace the domain")
    K = rtm.n_spectral_bins
    length(G_all) == K || error("G has $(length(G_all)) spectral bins but the domain has $K; retrace")
    Fv = Vector{AbstractMatrix}(F_matrices)

    ns = length(rtm.surface_mapping)
    N  = ns + length(rtm.volume_mapping)
    A  = n_angular_bins(rtm.directional_model)
    n  = N * A

    boundary, temperatures, emissive = setupBoundaryConditions(rtm, Fv)

    # ---- per-bin blocks: factorised once -------------------------------------
    verbose && println("  $N elements × $A direction bins = $n unknowns per bin, $K spectral bins: factorising...")
    ops = Vector{DirectionalSystemOp}(undef, K)
    lus = Vector{Any}(undef, K)
    Ms  = Vector{Matrix{Float64}}(undef, K)
    Ws  = Vector{Matrix{Float64}}(undef, K)          # (𝐃ᵀ𝐃)⁻¹ 𝐌ᵀ
    Whs = Vector{Vector{Float64}}(undef, K)          # its product with the boundary vector
    bs  = Vector{Vector{Float64}}(undef, K)
    inner = Matrix{Float64}(I, N, N)
    for k in 1:K
        op, _, _, b, known_q = _directional_setup(rtm, G_all[k]; spectral_bin = k)
        coeff = ifelse.(known_q, 1.0, b)
        D, Mk = _directional_dense_blocks(op, coeff)
        Fk = lu!(D)                                  # 𝐃 is not needed after this
        Wk = Fk \ (Fk' \ Matrix(Mk'))
        inner .+= Mk * Wk
        ops[k] = op; lus[k] = Fk; Ms[k] = Mk; Ws[k] = Wk; Whs[k] = Wk * boundary; bs[k] = b
    end
    inner_factor = cholesky(Symmetric(inner))

    sol      = [zeros(n) for _ in 1:K]
    previous = [zeros(n) for _ in 1:K]
    u_list   = [zeros(n) for _ in 1:K]
    Mu_total = zeros(N)
    emitFrac = getWeightedEmissionFractions(rtm, temperatures)

    # total emissive power from the current solution: Σ_k (j_k − ρ_k ∘ g_k)
    function emission_from_solution()
        em = zeros(N)
        for k in 1:K
            X = reshape(sol[k], N, A)
            _incident!(ops[k], X)
            em .+= vec(sum(X, dims = 2)) .- bs[k] .* vec(sum(ops[k].g, dims = 2))
        end
        return max.(em, 10 * eps(Float64))
    end

    # one sweep: Woodbury solve for the given emission state, then the updated state
    function sweep!(emissive_in, emitFrac_in, temperatures_in, iter)
        fill!(Mu_total, 0.0)
        for k in 1:K
            e_k = emissive_in .* view(emitFrac_in, :, k)
            u_list[k] .= Whs[k] .+ (lus[k] \ _spread_emission(ops[k], e_k))
            mul!(Mu_total, Ms[k], u_list[k], 1.0, 1.0)
        end
        w = inner_factor \ Mu_total
        for k in 1:K
            sol[k] .= u_list[k] .- Ws[k] * w
        end
        em  = iter > 1 ? emission_from_solution() :
              updateSpectralEmission!(rtm, 1, Fv, zeros(K * N), emitFrac_in, temperatures_in, emissive_in)
        Tn  = updateTemperaturesSpectral!(rtm, em, getBinsEmissionFractions(rtm, temperatures_in))
        efn = getWeightedEmissionFractions(rtm, Tn)
        return em, efn, Tn
    end

    verbose && println("  Starting spectral-directional sweeps...")
    prev_error = Inf                  # convergence error of the previous sweep (floor detection)
    went_up    = falses(10)           # for each of the last 10 sweeps: did the convergence error increase?
    max_dT     = Inf                  # largest temperature change of the latest sweep [K]
    for iter in 1:max_iters
        T_in_sweep = copy(temperatures)                 # temperatures fed into this sweep, for the reported change
        emissive, emitFrac, temperatures = sweep!(emissive, emitFrac, temperatures, iter)
        max_dT = maximum(abs.(temperatures .- T_in_sweep))   # largest temperature change of this sweep [K]

        num = 0.0; den = 0.0
        for k in 1:K
            num += sum(abs2, sol[k] .- previous[k]); den += sum(abs2, sol[k])
            previous[k] .= sol[k]
        end
        convergence_error = sqrt(num / den)
        (verbose && iter % 20 == 0) && println("  Iteration $iter: convergence error = $convergence_error, largest temperature change = $max_dT K")

        # Floor detection
        went_up[mod1(iter, 10)] = convergence_error > prev_error
        prev_error = convergence_error
        at_floor = convergence_error < 4 * eps(Float64) || (convergence_error < 1e-10 && count(went_up) >= 4)

        done_tol = iter > 1 && convergence_error < convergence_tol
        done     = done_tol || (iter > 1 && at_floor)
        if done || iter == max_iters
            done || @warn "Warning: Maximum iterations reached, final errors: convergence error = $convergence_error, largest temperature change = $max_dT K"
            emissive, emitFrac, temperatures = sweep!(emissive, emitFrac, temperatures, max(iter, 2))
            emissive     = emission_from_solution()
            temperatures = updateTemperaturesSpectral!(rtm, emissive, getBinsEmissionFractions(rtm, temperatures))
            if done
                if done_tol
                    verbose && println("  Converged after $iter iterations: convergence error = $convergence_error, ",
                            "largest temperature change = $max_dT K")
                else
                    verbose && println("  Converged to rounding level after $iter iterations: convergence error = $convergence_error, ",
                            "largest temperature change = $max_dT K (the requested tolerance $convergence_tol is below what the arithmetic resolves)")
                end
            end
            break
        end
    end

    # ---- results --------------------------------------------------------------
    verbose && println("  Writing spectral results to mesh...")
    energy_error = zeros(K)
    J_out = Vector{Matrix{Float64}}(undef, K)
    for k in 1:K
        X = reshape(sol[k], N, A)
        _incident!(ops[k], X)
        j_bin   = vec(sum(X, dims = 2))
        g_bin   = vec(sum(ops[k].g, dims = 2))
        r_bin   = bs[k] .* g_bin
        g_a_bin = g_bin .- r_bin
        writeResultsToDomain!(rtm, j_bin, g_a_bin, r_bin; spectral_bin = k)
        energy_error[k] = sum(j_bin .- g_bin) / (sum(j_bin) > 1000 * eps(Float64) ? sum(j_bin) : 1.0)
        J_out[k] = copy(X)
    end
    writeTemperaturesHeatSources!(rtm, temperatures)
    rtm.J = J_out
    rtm.energy_error  = energy_error

    if verbose
        show(stdout, MIME"text/plain"(), rtm)
        println()
    end
    return nothing
end