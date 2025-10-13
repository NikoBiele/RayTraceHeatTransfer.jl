function (rtm::RayTracingMeshOptim{VPF,VVPF,MT,VT,DIII,DII,GRID})(rays_tot::P; method::Symbol=:exchange,
                                nudge=nothing, rtol=nothing) where {VPF,VVPF,MT,VT,DIII,DII,P<:Integer,GRID}
    
    # Extract floating point type from the mesh vertices (Point2{G} where G is the float type)
    println("element type of rtm.fine_mesh[1][1].T_in_g is $(eltype(rtm.fine_mesh[1][1].T_in_g))")
    FloatType = eltype(rtm.fine_mesh[1][1].T_in_g) # Gets G from PolyFace2D{G}
    
    # Set defaults based on the mesh's floating point precision
    if FloatType <: Measurement{G} where G<:AbstractFloat
        actual_nudge = nudge === nothing ? 100 * (eps(FloatType) ± eps(FloatType)) : nudge
        actual_rtol = rtol === nothing ? 10 * (eps(FloatType) ± eps(FloatType)) : rtol
    else
        actual_nudge = nudge === nothing ? 100 * eps(FloatType) : nudge
        actual_rtol = rtol === nothing ? 10 * eps(FloatType) : rtol
    end

    if method == :exchange
        exchange_ray_tracing!(rtm, rays_tot, actual_rtol, actual_nudge)
    elseif method == :direct
        direct_ray_tracing!(rtm, rays_tot, actual_nudge)
    else
        error("Unknown ray tracing method: $method, must be :exchange or :direct")
    end
end