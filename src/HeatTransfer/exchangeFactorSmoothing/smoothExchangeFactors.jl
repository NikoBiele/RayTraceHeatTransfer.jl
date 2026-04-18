# Step 1: F = E * F (row scaling by diagonal E)
function step1_scale_rows_inplace!(F::Matrix{T}, scaling_factors::Vector{T}) where T
    n = size(F, 1)
    @inbounds for i in 1:n
        scale = scaling_factors[i]
        for j in 1:n
            F[i, j] *= scale
        end
    end
end

# Step 1: F = E * F (row scaling), fused with convergence measure
# Computes d = ||X1 - X1'||_F where X1 = E*F, using skew-symmetry.
function step1_scale_rows_with_conv!(F::Matrix{T}, scaling_factors::Vector{T}) where T
    n = size(F, 1)
    d_squared = zero(T)
    @inbounds for i in 1:n
        scale_i = scaling_factors[i]
        for j in (i+1):n
            x_ij = scale_i * F[i, j]
            x_ji = scaling_factors[j] * F[j, i]
            diff = x_ij - x_ji
            d_squared += diff * diff
            F[i, j] = x_ij
            F[j, i] = x_ji
        end
        F[i, i] = scale_i * F[i, i]
    end
    return sqrt(2 * d_squared)
end

# Step 2: F = 0.5 * (F + F^T) - symmetrization IN-PLACE
function step2_symmetrize_inplace!(F::Matrix{T}) where T
    n = size(F, 1)
    @inbounds for i in 1:n
        for j in (i+1):n  # Only upper triangle (excluding diagonal)
            # Compute symmetric value
            val = 0.5 * (F[i, j] + F[j, i])
            # Assign to both positions
            F[i, j] = val
            F[j, i] = val
        end
        # Diagonal remains unchanged: F[i,i] = 0.5*(F[i,i] + F[i,i]) = F[i,i]
    end
end

# Step 3: F = E^(-1) * F (row scaling by inverse diagonal)
function step3_scale_rows_inverse_inplace!(F::Matrix{T}, scaling_factors::Vector{T}) where T
    n = size(F, 1)
    @inbounds for i in 1:n
        inv_scale = 1.0 / scaling_factors[i]
        for j in 1:n
            F[i, j] *= inv_scale
        end
    end
end

# Step 4: F = diag(F*1)^(-1) * F - row normalization
function step4_normalize_rows_inplace!(F::Matrix{T}) where T
    n = size(F, 1)
    @inbounds for i in 1:n
        # Compute row sum
        row_sum = zero(T)
        for j in 1:n
            row_sum += F[i, j]
        end
        
        # Normalize row (if sum > 0)
        if row_sum > 0
            inv_sum = 1.0 / row_sum
            for j in 1:n
                F[i, j] *= inv_sum
            end
        end
    end
end

# Parallel versions
function step1_scale_rows_inplace_parallel!(F::Matrix{T}, scaling_factors::Vector{T}) where T
    n = size(F, 1)
    @threads for i in 1:n
        scale = scaling_factors[i]
        @inbounds for j in 1:n
            F[i, j] *= scale
        end
    end
end

# Parallel version
function step1_scale_rows_with_conv_parallel!(F::Matrix{T}, scaling_factors::Vector{T}) where T
    n = size(F, 1)
    d_squared = Threads.Atomic{T}(zero(T))
    Threads.@threads for i in 1:n
        scale_i = scaling_factors[i]
        local_sum = zero(T)
        @inbounds for j in (i+1):n
            x_ij = scale_i * F[i, j]
            x_ji = scaling_factors[j] * F[j, i]
            diff = x_ij - x_ji
            local_sum += diff * diff
            F[i, j] = x_ij
            F[j, i] = x_ji
        end
        @inbounds F[i, i] = scale_i * F[i, i]
        Threads.atomic_add!(d_squared, local_sum)
    end
    return sqrt(2 * d_squared[])
end    

function step2_symmetrize_inplace_parallel!(F::Matrix{T}) where T
    n = size(F, 1)
    @threads for i in 1:n
        @inbounds for j in (i+1):n
            val = 0.5 * (F[i, j] + F[j, i])
            F[i, j] = val
            F[j, i] = val
        end
    end
end

function step3_scale_rows_inverse_inplace_parallel!(F::Matrix{T}, scaling_factors::Vector{T}) where T
    n = size(F, 1)
    @threads for i in 1:n
        inv_scale = 1.0 / scaling_factors[i]
        @inbounds for j in 1:n
            F[i, j] *= inv_scale
        end
    end
end

function step4_normalize_rows_inplace_parallel!(F::Matrix{T}) where T
    n = size(F, 1)
    @threads for i in 1:n
        row_sum = zero(T)
        @inbounds for j in 1:n
            row_sum += F[i, j]
        end
        
        if row_sum > 0
            inv_sum = 1.0 / row_sum
            @inbounds for j in 1:n
                F[i, j] *= inv_sum
            end
        end
    end
end

# Updated smoothing algorithm - reads extinction directly from faces for specific spectral bin
function smoothExchangeFactors!(F::Matrix{T}, rtm::RayTracingDomain2D, 
                                 rays_per_emitter::Int,
                                 spectral_bin::Int=1; 
                                 max_iterations::Int=10_000,
                                 tolerance=nothing,
                                 check_interval::Int=2,
                                 stagnation_threshold=1e-4,
                                 verbose::Bool=true) where {T}
    
    num_surfaces = length(rtm.surface_mapping)
    num_volumes = length(rtm.volume_mapping)
    
    # Extract surface areas, volumes, and local extinction for specific spectral bin
    surface_areas = Vector{T}()
    volumes = Vector{T}() 
    volume_betas = Vector{T}()  # Local extinction coefficients for this bin
    
    for (coarse_index, coarse_face) in enumerate(rtm.coarse_mesh)
        for (fine_index, fine_face) in enumerate(rtm.fine_mesh[coarse_index])
            for (wall_index, is_solid) in enumerate(fine_face.solidWalls)
                if is_solid
                    push!(surface_areas, fine_face.area[wall_index])
                end
            end
            push!(volumes, fine_face.volume)
            
            # Extract local extinction for specific spectral bin
            if isa(fine_face.kappa_g, Vector)
                local_beta = fine_face.kappa_g[spectral_bin] + fine_face.sigma_s_g[spectral_bin]
            else
                local_beta = fine_face.kappa_g + fine_face.sigma_s_g
            end
            push!(volume_betas, local_beta)
        end
    end
    
    # Call ultimate smoothing with bin-specific extinction
    return smoothExchangeFactorsUltimate!(F, surface_areas, volumes, volume_betas, 
                                           rays_per_emitter; 
                                           max_iterations=max_iterations, 
                                           tolerance=tolerance,
                                           check_interval=check_interval,
                                           stagnation_threshold=stagnation_threshold,
                                           verbose=verbose)
end

function smoothExchangeFactorsUltimate!(F::Matrix{T}, surface_areas::Vector{P}, 
                                         volumes::Vector{P}, volume_betas::Vector{P},
                                         rays_per_emitter::Int; 
                                         max_iterations::Int=10_000,
                                         tolerance=nothing,
                                         check_interval::Int=2,
                                         stagnation_threshold=1e-4,
                                         verbose::Bool=true) where {T, P}
    
    num_surfaces = length(surface_areas)
    num_volumes = length(volumes)
    num_elements = num_surfaces + num_volumes
    
    # Determine underlying numeric types
    NumericT = if T <: Measurement
        T.parameters[1]
    else
        T
    end
    
    NumericP = if P <: Measurement
        P.parameters[1]
    else
        P
    end
    
    # Extract numeric values from input matrix
    F_work = Matrix{NumericT}(undef, size(F))
    for i in 1:size(F, 1), j in 1:size(F, 2)
        F_work[i, j] = if T <: Measurement
            Measurements.value(F[i, j])
        else
            F[i, j]
        end
    end

    # Pre-compute scaling factors with local extinction
    scaling_factors = Vector{NumericT}(undef, num_elements)
    for i in 1:num_elements
        if i <= num_surfaces
            area_val = surface_areas[i]
            scaling_factors[i] = if area_val isa Measurement
                Measurements.value(area_val)
            else
                area_val
            end
        else
            vol_idx = i - num_surfaces
            vol_val = volumes[vol_idx]
            beta_val = volume_betas[vol_idx]  # Use local extinction
            
            vol_numeric = if vol_val isa Measurement
                Measurements.value(vol_val)
            else
                vol_val
            end
            beta_numeric = if beta_val isa Measurement
                Measurements.value(beta_val)
            else
                beta_val
            end
            scaling_factors[i] = 4 * beta_numeric * vol_numeric
        end
    end

    # Adaptive tolerance based on input type
    numeric_tolerance = if tolerance === nothing
        κ = maximum(scaling_factors) / minimum(scaling_factors)
        sqrt(eps(NumericT)) * κ * sqrt(num_elements / rays_per_emitter)
    else
        if typeof(tolerance) <: Measurement
            Measurements.value(tolerance)
        else
            NumericT(tolerance)
        end
    end

    # Choose parallel strategy
    use_parallel = num_elements > 1000 && Threads.nthreads() > 1
    
    verbose && println("Matrix size: $(num_elements)×$(num_elements)")
    verbose && println("Strategy: $(use_parallel ? "Parallel" : "Serial")")
    verbose && println("Tolerance: $numeric_tolerance")

    # Pre-loop convergence check (sufficient condition from Algorithm 1)
    if num_volumes > 0
        Es_max = maximum(view(scaling_factors, 1:num_surfaces))
        Eg_min = minimum(view(scaling_factors, num_surfaces+1:num_elements))
        check_passed = Es_max < Eg_min
        if !check_passed && verbose
            @warn "Algorithm 1 convergence check failed: max surface E ($Es_max) ≥ min gas E ($Eg_min); convergence not guaranteed, consider refining the mesh"
        end
    else
        E_max = maximum(scaling_factors)
        E_sum = sum(scaling_factors)
        check_passed = E_max < 0.5 * E_sum
        if !check_passed && verbose
            @warn "Algorithm 1 convergence check failed: max surface E ($E_max) ≥ half of total E ($(E_sum/2)); convergence not guaranteed, consider refining the mesh"
        end
    end
    
    # Dispatch (now both Step 1 variants)
    if use_parallel
        step1_plain! = step1_scale_rows_inplace_parallel!
        step1_conv!  = step1_scale_rows_with_conv_parallel!
        symmetrize!  = step2_symmetrize_inplace_parallel!
        scale_rows_inverse! = step3_scale_rows_inverse_inplace_parallel!
        normalize_rows!     = step4_normalize_rows_inplace_parallel!
    else
        step1_plain! = step1_scale_rows_inplace!
        step1_conv!  = step1_scale_rows_with_conv!
        symmetrize!  = step2_symmetrize_inplace!
        scale_rows_inverse! = step3_scale_rows_inverse_inplace!
        normalize_rows!     = step4_normalize_rows_inplace!
    end

    d = Inf
    d_prev = Inf
    k = 0
    stagnated = false

    while d > numeric_tolerance && k < max_iterations
        if k % check_interval == 0
            d_prev = d
            d = step1_conv!(F_work, scaling_factors)
        else
            step1_plain!(F_work, scaling_factors)
        end
        symmetrize!(F_work)
        scale_rows_inverse!(F_work, scaling_factors)
        normalize_rows!(F_work)

        # Stagnation check: after at least two d samples, if d isn't shrinking
        if k > 2 * check_interval && d > numeric_tolerance && isfinite(d_prev) && d_prev > 0
            if (d_prev - d) / d_prev < stagnation_threshold
                stagnated = true
                @warn "No convergence progress detected at iteration $k; norm(E*F-F'*E) = $d. Stopping."
                break
            end
        end
        k += 1
        verbose && (k == 1 || k % 10 == 0) && println("Iteration $k: norm(E*F-F'*E) = $d")
    end

    if !stagnated && d < numeric_tolerance
        verbose && println("Converged after $k iterations. norm(E*F-F'*E) = $d")
    elseif !stagnated && k >= max_iterations && d > numeric_tolerance
        @warn "Max iterations ($max_iterations) reached without convergence. norm(E*F-F'*E) = $d"
    end

    # Handle uncertainty calculation if needed
    if T <: Measurement
        F_smooth_uncertain = Matrix{Measurement{NumericT}}(undef, size(F_work))
        @inbounds for i in 1:size(F_work, 1)
            for j in 1:size(F_work, 2)
                uncertainty = sqrt(F_work[i, j] / rays_per_emitter)
                F_smooth_uncertain[i, j] = F_work[i, j] ± uncertainty
            end
        end

        return F_smooth_uncertain
    end

    return F_work
end