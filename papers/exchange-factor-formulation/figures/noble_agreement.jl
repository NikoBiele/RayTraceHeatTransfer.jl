# Numerical agreement with Noble's matrix formulation of Hottel's zonal method (no figure). A unit square with 21 × 21
# cells, β = 1 m⁻¹, black walls with the bottom wall at 1000 K and the others at 0 K, and the gas in radiative
# equilibrium. For each single scattering albedo, the proposed formulation is solved, the blackbody emissive powers are
# recovered from its solution, and Noble's net source terms Q (surfaces) and S (gas) are computed from them. The
# largest difference from the source terms of the proposed formulation is printed, relative to the emitted power.
# The equivalence requires reciprocity, which ray-traced exchange factors satisfy only up to their noise, so the
# smoothed exchange factors (exactly reciprocal) are used here.

using RayTraceHeatTransfer                         # the GERT solver
using GeometryBasics, StaticArrays                 # points and fixed-size vectors for the geometry
using LinearAlgebra                                # Diagonal, inv, transpose

T_hot  = 1000.0                                    # temperature of the hot wall [K]
beta   = 1.0                                       # extinction coefficient [m⁻¹]
N      = 21                                        # cells per side
albedos = 0.0:0.1:0.9                              # single scattering albedos; at ω = 1 the gas has no blackbody emissive power in Noble's variables

vertices   = SVector(Point2(0.0, 0.0), Point2(1.0, 0.0), Point2(1.0, 1.0), Point2(0.0, 1.0))   # unit square, walls ordered bottom, right, top, left
solidWalls = SVector(true, true, true, true)       # all four walls are solid
face = PolyVolume2D{Float64}(vertices, solidWalls, 1, beta, 0.0)   # β = 1; the absorption–scattering split is set per albedo below
face.T_in_w  = [T_hot, 0.0, 0.0, 0.0]              # bottom wall hot, the others at 0 K
face.epsilon = [1.0, 1.0, 1.0, 1.0]                # black walls
face.T_in_g  = -1.0                                # gas temperature unknown ...
face.q_in_g  = 0.0                                 # ... in radiative equilibrium

mesh = RayTracingDomain2D([face], [(N, N)]; verbose = false)        # N × N cells
Ns = length(mesh.surface_mapping)                  # number of surface elements, numbered first
Ng = length(mesh.volume_mapping)                   # number of gas elements, numbered after the surfaces
mesh((Ns + Ng) * 2^16; method = :exchange, sampler = :sobol, verbose = false)   # exchange factors by ray tracing, 2¹⁶ rays per element
smooth!(mesh; verbose = false)                     # exactly reciprocal exchange factors, as the equivalence requires

A = Float64.(mesh.surface_areas)                   # surface areas (wall lengths in 2D) [m]
V = Float64.(mesh.volumes)                         # gas volumes (cell areas in 2D) [m²]
w = vcat(A, 4 .* beta .* V)                        # emission capacities: A for surfaces, 4βV for gas
s = 1:Ns                                           # indices of the surface elements
g = Ns+1:Ns+Ng                                     # indices of the gas elements

# Largest reciprocity defect of the exchange factors, relative to the largest direct exchange area
reciprocity_defect(F) = (WF = Diagonal(w) * Matrix(F); maximum(abs.(WF .- transpose(WF))) / maximum(abs.(WF)))
println("relative reciprocity defect: raw F ", reciprocity_defect(mesh.F_raw), ", smoothed F ", reciprocity_defect(mesh.F_smooth))

F  = Matrix(mesh.F_smooth)                         # dense exchange factors
WF = Diagonal(w) * F                               # direct exchange areas
ss, sg = WF[s, s], WF[s, g]                        # Noble's blocks: surface–surface and surface–gas ...
gs, gg = WF[g, s], WF[g, g]                        # ... gas–surface and gas–gas
ε = ones(Ns)                                       # surface emissivities (black walls)
ρ = 1 .- ε                                         # surface reflectivities

for omega in albedos
    for cell in mesh.fine_mesh[1]                  # same extinction, so the same exchange factors ...
        cell.kappa_g   = (1 - omega) * beta        # ... with absorption κ = (1 − ω)β ...
        cell.sigma_s_g = omega * beta              # ... and scattering σₛ = ωβ
    end
    solveEquilibrium!(mesh, mesh.F_smooth; verbose = false)   # the proposed formulation: one linear solve

    # Source terms of the proposed formulation from its solution j (power leaving each element)
    j = mesh.J                                     # radiant powers [W]
    b = vcat(ρ, fill(omega, Ng))                   # reflectivity of surfaces, albedo of gas
    incident = transpose(F) * j                    # power arriving at each element [W]
    e = j .- b .* incident                         # emitted power [W]
    q = e .- (1 .- b) .* incident                  # net source: emitted minus absorbed [W]

    # Blackbody emissive powers, recovered as e / ((1 − b) w): εA for surfaces, 4κV for gas [W/m²]
    E  = e[s] ./ (ε .* A)
    Eg = e[g] ./ ((1 - omega) .* 4 .* beta .* V)

    # Noble's intermediate matrices and total exchange areas. Noble writes P and R as inverses; here each operator is
    # factorized once and applied by linear solves, so P⁻¹X becomes P_lu \ X and R⁻¹X becomes R_lu \ X.
    P_lu = lu(Diagonal(4 .* beta .* V) .- omega .* gg)                      # gas multiple-scattering operator, factorized
    L    = sg * (P_lu \ gs)                                                 # surface–gas–surface exchange via scattering
    R_lu = lu(Diagonal(A) .- (ss .+ omega .* L) * Diagonal(ρ))              # surface multiple-reflection operator, factorized
    K    = sg * (P_lu \ Matrix(Diagonal(4 .* beta .* V)))                   # surface–gas coupling
    SS   = Diagonal(ε .* A) * (R_lu \ (ss .+ omega .* L)) * Diagonal(ε)     # surface–surface total exchange areas
    SG   = (1 - omega) .* Diagonal(ε .* A) * (R_lu \ K)                     # surface–gas total exchange areas
    GS   = transpose(SG)                                                    # gas–surface, by reciprocity
    GG   = (1 - omega)^2 .* (Diagonal(4 .* beta .* V) * (P_lu \ gg) .+ transpose(K) * Diagonal(ρ) * (R_lu \ K))   # gas–gas total exchange areas
    
    # Noble's net source terms: Q emitted minus absorbed (surfaces), S absorbed minus emitted (gas)
    Q = ε .* A .* E .- SS * E .- SG * Eg
    S = GG * Eg .+ GS * E .- (1 - omega) .* 4 .* beta .* V .* Eg

    scale = sum(e[s])                              # power emitted by the walls, for a relative measure
    println("ω = $omega: largest relative difference, surfaces ", maximum(abs.(Q .- q[s])) / scale,
            ", gas ", maximum(abs.(S .+ q[g])) / scale)   # S has the opposite sign convention to q
end

## Hottel's zonal balances with the raw, non-reciprocal exchange factors (Theorem 5.1, Remark 5.2).
# Hottel and Sarofim define the direct exchange area from the emitter: s_i s_k = w_i F_ik is the power leaving i
# that arrives at k, so the matrix X with X_ik = w_i F_ik is summed over the emitters i to give the incident power
# on k. Hottel's balances are solved directly for the leaving flux densities with this X, and the net source terms
# are compared with those of the proposed formulation and summed. For reference, the same is done with the
# transposed matrix, which the symmetric notation would permit only for reciprocal F.

omega = 0.5                                        # one albedo suffices: with black walls and the gas in equilibrium, M does not depend on it
for cell in mesh.fine_mesh[1]                      # set the absorption–scattering split for this albedo
    cell.kappa_g   = (1 - omega) * beta            # absorption κ = (1 − ω)β
    cell.sigma_s_g = omega * beta                  # scattering σₛ = ωβ
end
solveEquilibrium!(mesh, mesh.F_raw; verbose = false)   # the proposed formulation with the raw exchange factors
F_raw = Matrix(mesh.F_raw)                         # dense raw exchange factors, as traced
j_raw = mesh.J                                     # radiant powers of the proposed formulation [W]
b     = vcat(ρ, fill(omega, Ng))                   # reflectivity of surfaces, albedo of gas
g_raw = transpose(F_raw) * j_raw                   # incident powers of the proposed formulation [W]
e_raw = j_raw .- b .* g_raw                        # emitted powers [W]
q_raw = e_raw .- (1 .- b) .* g_raw                 # net source terms of the proposed formulation [W]
Eb    = e_raw ./ ((1 .- b) .* w)                   # blackbody emissive powers, e / (εA) for surfaces and e / (4κV) for gas [W/m²]

# Hottel's balances with direct exchange areas X, solved for the leaving flux densities u:
# u = (1 − b) E_b + b H with incident flux densities H = X u / w; net source q = w (u − H) = (1 − b) w (E_b − H)
function hottel_balances(X, Eb, b, w)
    T = Diagonal(1 ./ w) * X                       # incident flux density per unit leaving flux density
    u = (I - Diagonal(b) * T) \ ((1 .- b) .* Eb)   # leaving flux densities [W/m²]
    H = T * u                                      # incident flux densities [W/m²]
    return w .* (u .- H)                           # net source terms [W]
end

X_hottel     = transpose(Diagonal(w) * F_raw)      # Hottel's definition: X[k, i] = w_i F_ik, so X * u sums over the emitters i
X_transposed = Diagonal(w) * F_raw                 # the transposed areas: X[k, i] = w_k F_ki, equal to the above only for reciprocal F
scale = sum(e_raw[s])                              # power emitted by the walls, for relative measures
println("proposed formulation, raw F: relative energy conservation error ", abs(sum(q_raw)) / scale)
for (name, X) in (("Hottel's definition", X_hottel), ("transposed areas", X_transposed))
    q_hottel = hottel_balances(X, Eb, b, w)        # net source terms of Hottel's balances [W]
    println("Hottel's balances, raw F, ", name, ": relative energy conservation error ", abs(sum(q_hottel)) / scale,
            ", largest relative difference from the proposed formulation ", maximum(abs.(q_hottel .- q_raw)) / scale)
end