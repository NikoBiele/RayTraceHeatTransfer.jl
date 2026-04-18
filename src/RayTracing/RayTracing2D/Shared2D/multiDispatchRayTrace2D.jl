function (rtm::RayTracingDomain2D{VPF,VVPF,MT,VT,DIII,DII,GRID})(rays_tot::P; method::Symbol=:exchange,
                                nudge=nothing, smooth_tol=nothing,
                                check_interval=2, stagnation_threshold=1e-4, verbose=true) where
                                {VPF,VVPF,MT,VT,DIII,DII,P<:Integer,GRID}
    
    # Extract floating point type from the mesh vertices (Point2{G} where G is the float type)
    FloatType = eltype(rtm.fine_mesh[1][1].T_in_g) # Gets G from PolyVolume2D{G}
    
    # Set defaults based on the mesh's floating point precision
    if FloatType <: Measurement{G} where G<:AbstractFloat
        trace_nudge = nudge === nothing ? 100 * (eps(FloatType) ± eps(FloatType)) : nudge
        smooth_rtol = smooth_tol === nothing ? nothing : smooth_tol # 10 * (eps(FloatType) ± eps(FloatType))
    else
        trace_nudge = nudge === nothing ? 100 * eps(FloatType) : nudge
        smooth_rtol = smooth_tol === nothing ? nothing : smooth_tol # 10 * eps(FloatType)
    end

    if method == :exchange
        exchangeRayTracing!(rtm, rays_tot, smooth_rtol, trace_nudge,
                            check_interval, stagnation_threshold, verbose)
    elseif method == :direct
        directRayTracing!(rtm, rays_tot, trace_nudge, verbose)
    else
        error("Unknown ray tracing method: $method, must be :exchange or :direct")
    end
end