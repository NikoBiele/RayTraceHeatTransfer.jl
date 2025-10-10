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

# Convergence check: compute d = ||E*F - F^T*E||_F using only scalars
function compute_convergence_metric(F::Matrix{T}, scaling_factors::Vector{T}) where T
    n = size(F, 1)
    d_squared = zero(T)
    
    @inbounds for i in 1:n
        for j in 1:n
            # (E*F)[i,j] = scaling_factors[i] * F[i,j]
            EF_ij = scaling_factors[i] * F[i, j]
            
            # (F^T*E)[i,j] = F[j,i] * scaling_factors[j]
            FTE_ij = F[j, i] * scaling_factors[j]
            
            # ||E*F - F^T*E||_F^2 = sum_{i,j} (EF_ij - FTE_ij)^2
            diff = EF_ij - FTE_ij
            d_squared += diff * diff
        end
    end
    
    return sqrt(d_squared)
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

function compute_convergence_metric_parallel(F::Matrix{T}, scaling_factors::Vector{T}) where T
    n = size(F, 1)
    
    # Parallel reduction for d_squared
    d_squared = Threads.Atomic{T}(zero(T))
    
    @threads for i in 1:n
        local_sum = zero(T)
        @inbounds for j in 1:n
            EF_ij = scaling_factors[i] * F[i, j]
            FTE_ij = F[j, i] * scaling_factors[j]
            diff = EF_ij - FTE_ij
            local_sum += diff * diff
        end
        # Atomic add to global sum
        Threads.atomic_add!(d_squared, local_sum)
    end
    
    return sqrt(d_squared[])
end

# Updated smoothing algorithm - reads extinction directly from faces for specific spectral bin
function smooth_exchange_factors!(F::Matrix{T}, rtm::RayTracingMeshOptim, 
                                 rays_per_emitter::Int,
                                 spectral_bin::Int=1; 
                                 max_iterations::Int=1000, tolerance=nothing) where {T}
    
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
    return smooth_exchange_factors_ultimate!(F, surface_areas, volumes, volume_betas, 
                                           rays_per_emitter; 
                                           max_iterations=max_iterations, 
                                           tolerance=tolerance)
end

function smooth_exchange_factors_ultimate!(F::Matrix{T}, surface_areas::Vector{P}, 
                                         volumes::Vector{P}, volume_betas::Vector{P},
                                         rays_per_emitter::Int; 
                                         max_iterations::Int=1000, tolerance=nothing) where {T, P}
    
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
    
    # Adaptive tolerance based on input type
    numeric_tolerance = if tolerance === nothing
        eps(NumericT)
    else
        if typeof(tolerance) <: Measurement
            Measurements.value(tolerance)
        else
            tolerance
        end
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
    
    # Choose parallel strategy
    use_parallel = num_elements > 1000 && Threads.nthreads() > 1
    
    println("Matrix size: $(num_elements)×$(num_elements)")
    println("Strategy: $(use_parallel ? "Parallel" : "Serial")")
    println("Tolerance: $numeric_tolerance")

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
    
    # Choose function implementations based on parallel strategy
    if use_parallel
        scale_rows! = step1_scale_rows_inplace_parallel!
        symmetrize! = step2_symmetrize_inplace_parallel!
        scale_rows_inverse! = step3_scale_rows_inverse_inplace_parallel!
        normalize_rows! = step4_normalize_rows_inplace_parallel!
        convergence_metric = compute_convergence_metric_parallel
    else
        scale_rows! = step1_scale_rows_inplace!
        symmetrize! = step2_symmetrize_inplace!
        scale_rows_inverse! = step3_scale_rows_inverse_inplace!
        normalize_rows! = step4_normalize_rows_inplace!
        convergence_metric = compute_convergence_metric
    end
    
    for iteration in 1:max_iterations
        # Step 1: F = E * F (overwrites F_work)
        scale_rows!(F_work, scaling_factors)
        
        # Step 2: F = 0.5 * (F + F^T) (overwrites F_work)
        symmetrize!(F_work)
        
        # Step 3: F = E^(-1) * F (overwrites F_work)
        scale_rows_inverse!(F_work, scaling_factors)
        
        # Step 4: F = row_normalize(F) (overwrites F_work)
        normalize_rows!(F_work)

        # Check convergence: d = ||E*F - F^T*E||_F
        d = convergence_metric(F_work, scaling_factors)

        if d < numeric_tolerance
            println("Converged after $iteration iterations. d = $d")
            break
        end

        # Progress reporting
        if iteration == 1 || iteration % 10 == 0
            println("Iteration $iteration: d = $d")
        end

        if iteration == max_iterations
            println("Warning: Maximum iterations reached. Final d = $d")
        end
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