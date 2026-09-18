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

# Cell-averaged (Galerkin) kernel on optical-depth edges τ_e (Nx+1).
function lbl_kernel!(K::Matrix{Float64}, w1::Vector{Float64}, w2::Vector{Float64},
                     τ_e::Vector{Float64})
    Nx = length(w1)
    τ_L = τ_e[end]
    @inbounds for i in 1:Nx
        a, b = τ_e[i], τ_e[i+1]
        h = b - a
        if h <= 0
            w1[i] = lbl_E2(a); w2[i] = lbl_E2(τ_L - a)
            for j in 1:Nx
                c, d = τ_e[j], τ_e[j+1]
                K[i, j] = j == i ? 0.0 :
                          (c >= b ? lbl_E2(c - b) - lbl_E2(d - b) : lbl_E2(a - d) - lbl_E2(a - c))
            end
            continue
        end
        w1[i] = (lbl_E3(a) - lbl_E3(b)) / h
        w2[i] = (lbl_E3(τ_L - b) - lbl_E3(τ_L - a)) / h
        for j in 1:Nx
            c, d = τ_e[j], τ_e[j+1]
            if j == i
                K[i, j] = 2 * (h - 0.5 + lbl_E3(h)) / h
            elseif c >= b
                K[i, j] = (lbl_E3(c - b) - lbl_E3(d - b) - lbl_E3(c - a) + lbl_E3(d - a)) / h
            else
                K[i, j] = (lbl_E3(a - d) - lbl_E3(b - d) - lbl_E3(a - c) + lbl_E3(b - c)) / h
            end
        end
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