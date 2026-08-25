"""
Pretty printing for the domain types.

Without these, echoing a domain in the REPL dumps the entire nested mesh and
both exchange factor matrices. These `show` methods give a compact summary
instead, and report which stage of the workflow the domain has reached:

    mesh  ->  exchange factors  ->  smooth  ->  solve

`show(io, ::MIME"text/plain", d)` is the multi-line REPL form.
`show(io, d)` is the one-line form used inside containers.

Neither method may throw: a domain that is malformed or half-built must still
print. Field access is guarded accordingly.
"""

using SparseArrays

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

# F_raw / F_smooth are allocated with `undef` contents at construction, so their
# presence says nothing about whether they hold results. Uninitialised memory is
# either NaN or denormal (~1e-310); a real exchange factor row sums to ~1. One
# row is enough to tell them apart, and it is O(N) rather than O(N^2).
function _computed(F::AbstractMatrix)
    (isempty(F) || size(F, 1) == 0) && return false
    s = sum(abs, @view F[1, :])
    return isfinite(s) && s > 1e-100
end
_computed(F::AbstractVector) = !isempty(F) && all(_computed, F)
_computed(::Nothing) = false
_computed(::Any) = false

function _stage(d)
    raw    = isdefined(d, :F_raw)    && _computed(getfield(d, :F_raw))
    smooth = isdefined(d, :F_smooth) && _computed(getfield(d, :F_smooth))
    raw    || return :meshed
    smooth || return :traced
    return :smoothed
end

_mdesc(F::AbstractMatrix) =
    string(size(F, 1), "×", size(F, 2), issparse(F) ? " sparse" : " dense")
_mdesc(F::AbstractVector) = isempty(F) ? "empty" :
    string(length(F), " × (", _mdesc(first(F)), ")")
_mdesc(::Any) = "—"

function _fmt_range(lo, hi)
    lo == hi && return string(round(lo; sigdigits = 4))
    return string(round(lo; sigdigits = 4), " … ", round(hi; sigdigits = 4))
end

_fmt_err(::Nothing) = "not solved"
_fmt_err(e::Real)   = string(round(e; sigdigits = 3))
_fmt_err(e::AbstractVector) = isempty(e) ? "not solved" :
    string(round(maximum(abs, e); sigdigits = 3), " (max over ", length(e), " bins)")
_fmt_err(::Any) = "—"

function _spectral_line(d)
    mode = getfield(d, :spectral_mode)
    mode === :grey && return "grey"
    n = getfield(d, :n_spectral_bins)
    s = string(mode, ", ", n, " bins")
    lims = getfield(d, :wavelength_band_limits)
    if lims isa AbstractVector && length(lims) >= 2
        s *= string(", λ ∈ [", round(first(lims); sigdigits = 3), ", ",
                    round(last(lims);  sigdigits = 3), "]")
    end
    return s
end

# number of distinct vertices in a face row: a repeated cyclic pair means the
# face was written as a degenerate quad and meshed as a triangle
function _arity(row)
    n = length(row)
    dup = count(i -> row[i] == row[mod1(i + 1, n)], 1:n)
    return n - dup
end

function _topology(faces::AbstractMatrix)
    nf = size(faces, 1)
    ar = [_arity(@view faces[f, :]) for f in 1:nf]
    tri  = count(==(3), ar)
    quad = count(==(4), ar)
    other = nf - tri - quad
    parts = String[]
    quad  > 0 && push!(parts, string(quad, " quad"))
    tri   > 0 && push!(parts, string(tri, " tri"))
    other > 0 && push!(parts, string(other, " other"))
    return isempty(parts) ? "" : string(" (", join(parts, ", "), ")")
end

# ---------------------------------------------------------------------------
# RayTracingDomain2D
# ---------------------------------------------------------------------------

function Base.show(io::IO, d::RayTracingDomain2D)
    ncoarse = length(d.coarse_mesh)
    nvol    = length(d.volume_mapping)
    nsurf   = length(d.surface_mapping)
    print(io, "RayTracingDomain2D(", ncoarse, " faces, ", nvol, " volumes, ",
              nsurf, " surfaces, ", d.spectral_mode, ")")
end

function Base.show(io::IO, ::MIME"text/plain", d::RayTracingDomain2D)
    get(io, :compact, false) && return show(io, d)

    ncoarse = length(d.coarse_mesh)
    nvol    = length(d.volume_mapping)
    nsurf   = length(d.surface_mapping)
    stage   = _stage(d)

    println(io, "RayTracingDomain2D")

    # geometry
    per = [length(f) for f in d.fine_mesh]
    geo = string(ncoarse, " coarse face", ncoarse == 1 ? "" : "s", " → ",
                 nvol, " volumes, ", nsurf, " surfaces")
    if !isempty(per) && minimum(per) != maximum(per)
        geo *= string("  (", minimum(per), "–", maximum(per), " per face)")
    end
    println(io, "  geometry   ", geo)
    d.surfaces_only && println(io, "             surfaces only (no participating medium)")

    # boundary conditions
    let npT = 0, nuT = 0
        for ((ci, fi, wi), _) in d.surface_mapping
            (d.fine_mesh[ci][fi].T_in_w[wi] < -0.1 ? (nuT += 1) : (npT += 1))
        end
        for ((ci, fi), _) in d.volume_mapping
            (d.fine_mesh[ci][fi].T_in_g < -0.1 ? (nuT += 1) : (npT += 1))
        end
        println(io, "  boundary   ", npT, " prescribed T, ", nuT,
                    " prescribed source")
    end

    println(io, "  spectral   ", _spectral_line(d))

    # workflow stage
    if stage === :meshed
        println(io, "  exchange   not computed — call domain(N_rays; method = :exchange)")
    else
        println(io, "  exchange   F_raw     ", _mdesc(d.F_raw))
        if stage === :traced
            println(io, "             F_smooth  not computed — call smooth!(domain)")
        else
            println(io, "             F_smooth  ", _mdesc(d.F_smooth))
        end
    end

    print(io, "  energy     ", _fmt_err(d.energy_error))
end

# ---------------------------------------------------------------------------
# ViewFactorDomain3D
# ---------------------------------------------------------------------------

function Base.show(io::IO, d::ViewFactorDomain3D)
    nf   = length(d.facesMesh)
    nsub = sum(length(f.subFaces) for f in d.facesMesh; init = 0)
    print(io, "ViewFactorDomain3D(", nf, " faces, ", nsub, " subfaces, ",
              d.spectral_mode, ")")
end

function Base.show(io::IO, ::MIME"text/plain", d::ViewFactorDomain3D)
    get(io, :compact, false) && return show(io, d)

    nf    = length(d.facesMesh)
    per   = [length(f.subFaces) for f in d.facesMesh]
    nsub  = sum(per; init = 0)
    stage = _stage(d)

    println(io, "ViewFactorDomain3D")

    # geometry — note the topology breakdown, since triangles and quads give
    # different subface counts at the same Ndim
    println(io, "  geometry   ", nf, " face", nf == 1 ? "" : "s",
                _topology(d.faces), " → ", nsub, " subfaces  (Ndim = ", d.Ndims, ")")
    if !isempty(per) && minimum(per) != maximum(per)
        println(io, "             ", minimum(per), "–", maximum(per),
                    " subfaces per face (ragged index map)")
    end

    # boundary conditions
    let npT = 0, nuT = 0, lo = Inf, hi = -Inf
        for f in d.facesMesh, sf in f.subFaces
            T = sf.T_in_w
            if T isa Real && T >= -0.1
                npT += 1; lo = min(lo, T); hi = max(hi, T)
            else
                nuT += 1
            end
        end
        s = string(npT, " prescribed T, ", nuT, " prescribed source")
        npT > 0 && (s *= string("  (T = ", _fmt_range(lo, hi), " K)"))
        println(io, "  boundary   ", s)
    end

    println(io, "  spectral   ", _spectral_line(d))

    if stage === :meshed
        println(io, "  exchange   not computed — call domain(; parallel = true)")
    else
        println(io, "  exchange   F_raw     ", _mdesc(d.F_raw))
        if stage === :traced
            println(io, "             F_smooth  not computed — call smooth!(domain)")
        else
            println(io, "             F_smooth  ", _mdesc(d.F_smooth))
        end
    end

    print(io, "  energy     ", _fmt_err(d.energy_error))
end

function Base.show(io::IO, f::PolyVolume2D{G}) where {G}
    print(io, "PolyVolume2D{", G, "}(", length(f.vertices), " vertices, ",
              count(f.solidWalls), "/", length(f.solidWalls), " solid, ",
              length(f.subVolumes), " subvolumes)")
end

function Base.show(io::IO, f::PolyFace3D{G}) where {G}
    nv = f.vertices === nothing ? 0 : length(f.vertices)
    ns = f.subFaces === nothing ? 0 : length(f.subFaces)
    print(io, "PolyFace3D{", G, "}(", nv, " vertices, ", ns, " subfaces)")
end