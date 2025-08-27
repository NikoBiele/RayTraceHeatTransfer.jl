function check_energy_conservation(F::Matrix{T}, tol::P=1e-12) where {T, P}
    num_emitters = size(F, 1)
    conservation_satisfied = true
    max_error = 0.0
    
    for i in 1:num_emitters
        row_sum = sum(F[i, :])
        error = abs(row_sum - 1.0)
        if error > max_error
            max_error = error
        end
        if !isapprox(row_sum, 1.0, atol=tol)
            conservation_satisfied = false
            println("Energy conservation violated for emitter $i: sum = $row_sum")
        end
    end
    
    if conservation_satisfied
        println("Energy conservation satisfied for all emitters within tolerance.")
    else
        println("Energy conservation not satisfied. Maximum error: $max_error")
    end
    
    return conservation_satisfied, max_error
end