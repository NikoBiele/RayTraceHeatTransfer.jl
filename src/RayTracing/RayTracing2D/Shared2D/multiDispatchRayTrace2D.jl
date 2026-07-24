function (rtm::RayTracingDomain2D{VPF,VVPF,MT,VT,DIII,DII,GRID})(rays_tot::P; method::Symbol=:exchange,
                                nudge=nothing, k_dykstra::Union{Nothing,Int}=nothing,
                                max_iters::P=1_000, verbose=true) where
                                {VPF,VVPF,MT,VT,DIII,DII,P<:Integer,GRID}
    
    # Extract floating point type from the mesh vertices (Point2{G} where G is the float type)
    FloatType = eltype(rtm.fine_mesh[1][1].T_in_g) # Gets G from PolyVolume2D{G}
    
    # Set defaults based on the mesh's floating point precision
    trace_nudge = nudge === nothing ? 10_000 * eps(FloatType) : nudge

    if method == :exchange
        exchangeRayTracing!(rtm, rays_tot, trace_nudge, max_iters, k_dykstra, verbose)
    elseif method == :direct
        directRayTracing!(rtm, rays_tot, trace_nudge, verbose)
    else
        error("Unknown ray tracing method: $method, must be :exchange or :direct")
    end
end