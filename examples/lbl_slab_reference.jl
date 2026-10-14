using RayTraceHeatTransfer
using RayTraceHeatTransfer: n_pieces
using GeometryBasics, StaticArrays
using SpecialFunctions      # expint
using LinearAlgebra
using Random
using Base.Threads

const LBL_σ  = 5.670374419e-8
const LBL_C1 = 1.191042972e-16
const LBL_C2 = 1.4387769e-2

"""
    lbl_slab_equilibrium(λ, κ, L, T1, T2; Nx, tol, max_iters)

Radiative-equilibrium cell temperatures of a slab of thickness `L` between
black plates at `T1`, `T2` for the spectrum `κ` sampled at `λ`.
Returns `(T, q, iters)`: cell temperatures, net flux from plate 1 to plate 2,
Newton iterations.
"""
function lbl_slab_equilibrium(λ::AbstractVector, κ::AbstractVector, L::Real, T1::Real, T2::Real;
                              Nx::Int = 40, tol::Real = 1e-8, max_iters::Int = 100)
    Nλ = length(λ)
    Nλ == length(κ) || throw(ArgumentError("λ and κ must have the same length"))

    wλ = zeros(Nλ)
    wλ[1] = (λ[2] - λ[1]) / 2
    wλ[end] = (λ[end] - λ[end-1]) / 2
    @inbounds for k in 2:Nλ-1
        wλ[k] = (λ[k+1] - λ[k-1]) / 2
    end

    Δx  = L / Nx
    S_e = collect(range(0.0, L, length = Nx + 1))

    T = fill((T1 + T2) / 2, Nx)
    A = zeros(Nx)
    J = zeros(Nx, Nx)

    nt = nthreads()
    chunks = [round(Int, (c - 1) * Nλ / nt) + 1 : round(Int, c * Nλ / nt) for c in 1:nt]
    Ks  = [zeros(Nx, Nx) for _ in 1:nt]
    Js  = [zeros(Nx, Nx) for _ in 1:nt]
    w1s = [zeros(Nx) for _ in 1:nt]
    w2s = [zeros(Nx) for _ in 1:nt]
    τes = [zeros(Nx + 1) for _ in 1:nt]
    As  = [zeros(Nx) for _ in 1:nt]
    Bs  = [zeros(Nx) for _ in 1:nt]
    dBs = [zeros(Nx) for _ in 1:nt]

    emission(Ti)  = 4π * sum(κ[k] * lbl_planck(λ[k], Ti) * wλ[k] for k in 1:Nλ)
    demission(Ti) = 4π * sum(κ[k] * lbl_dplanck(λ[k], Ti) * wλ[k] for k in 1:Nλ)

    r  = zeros(Nx)
    dE = zeros(Nx)
    iters = 0
    for it in 1:max_iters
        iters = it
        foreach(a -> fill!(a, 0.0), As)
        foreach(j -> fill!(j, 0.0), Js)

        @threads for c in 1:nt
            K, Jℓ, w1, w2, τ_e = Ks[c], Js[c], w1s[c], w2s[c], τes[c]
            Aℓ, B, dB = As[c], Bs[c], dBs[c]
            for k in chunks[c]
                κk = κ[k]
                τ_e .= κk .* S_e
                lbl_kernel!(K, w1, w2, τ_e)
                B1, B2 = lbl_planck(λ[k], T1), lbl_planck(λ[k], T2)
                @inbounds for j in 1:Nx
                    B[j]  = lbl_planck(λ[k], T[j])
                    dB[j] = lbl_dplanck(λ[k], T[j])
                end
                @inbounds for i in 1:Nx
                    Gi  = B1 * w1[i] + B2 * w2[i]
                    c_i = 2π * κk * wλ[k]
                    for j in 1:Nx
                        Gi += B[j] * K[i, j]
                        Jℓ[i, j] += c_i * K[i, j] * dB[j]
                    end
                    Aℓ[i] += c_i * Gi
                end
            end
        end
        fill!(A, 0.0); fill!(J, 0.0)
        for c in 1:nt
            A .+= As[c]
            J .+= Js[c]
        end

        @threads for i in 1:Nx
            r[i]  = emission(T[i]) - A[i]
            dE[i] = demission(T[i])
        end
        δ = -((Diagonal(dE) - J) \ r)
        α = minimum(δ[i] < 0 ? min(1.0, 0.9 * T[i] / -δ[i]) : 1.0 for i in 1:Nx)
        T .+= α .* δ
        maximum(abs.(α .* δ)) < tol && break
    end

    qλ = zeros(Nλ)
    @threads for k in 1:Nλ
        κk  = κ[k]
        τ_e = κk .* S_e
        acc = lbl_planck(λ[k], T2) * lbl_E3(τ_e[end])
        for j in 1:Nx
            acc += lbl_planck(λ[k], T[j]) * (lbl_E3(τ_e[j]) - lbl_E3(τ_e[j+1]))
        end
        qλ[k] = π * (lbl_planck(λ[k], T1) - 2 * acc) * wλ[k]
    end
    return T, sum(qλ), iters
end

# Cell-averaged (Galerkin) kernel on uniformly spaced optical-depth edges τ_e (Nx+1).
# On uniform cells the kernel depends only on |i − j| (Toeplitz), and every entry is
# a difference of E₃ at the edges, so Nx+1 evaluations of E₃ build all of it.
function lbl_kernel!(K::Matrix{Float64}, w1::Vector{Float64}, w2::Vector{Float64},
                     τ_e::Vector{Float64})
    Nx = length(w1)
    h  = τ_e[2] - τ_e[1]                         # optical thickness of one cell
    if h <= 0                                    # transparent: no exchange with the gas
        fill!(K, 0.0)                            # no gas-to-gas coupling
        fill!(w1, 1.0)                           # E₂(0) = 1: plate 1 seen unattenuated
        fill!(w2, 1.0)                           # E₂(0) = 1: plate 2 seen unattenuated
        return nothing
    end
    e = [lbl_E3(τ) for τ in τ_e]                 # E₃ at every edge, e[k] = E₃((k−1)h)
    @inbounds for i in 1:Nx
        w1[i] = (e[i] - e[i+1]) / h              # cell i seen from plate 1 (at τ = 0)
        w2[i] = (e[Nx+1-i] - e[Nx+2-i]) / h      # cell i seen from plate 2 (at τ = τ_L)
    end
    g = zeros(Nx)                                # kernel as a function of m = |i − j|
    g[1] = 2 * (h - 0.5 + e[2]) / h              # self-exchange of a cell (m = 0)
    @inbounds for m in 1:Nx-1
        g[m+1] = (e[m] - e[m+1] - e[m+1] + e[m+2]) / h   # cells m apart
    end
    @inbounds for j in 1:Nx, i in 1:Nx
        K[i, j] = g[abs(i - j) + 1]              # fill the Toeplitz matrix
    end
    return nothing
end

@inline function lbl_E2(u::Float64)
    u <= 0 && return 1.0
    u > 700 && return 0.0
    return exp(-u) - u * expint(u)
end
@inline lbl_E3(u::Float64) = u <= 0 ? 0.5 : (u > 700 ? 0.0 : (exp(-u) - u * lbl_E2(u)) / 2)

lbl_planck(λ, T)  = LBL_C1 / (λ^5 * expm1(LBL_C2 / (λ * T)))
lbl_dplanck(λ, T) = (x = LBL_C2 / (λ * T); lbl_planck(λ, T) * (x / T) / (-expm1(-x)))