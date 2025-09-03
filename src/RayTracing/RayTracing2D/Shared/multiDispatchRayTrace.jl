function (rtm::RayTracingMeshOptim{VPF,VVPF,MT,VT,DIII,DII,BOO,FLOA,GRID})(rays_tot::P; method::Symbol=:exchange, nudge=nothing, rtol=nothing, uncertain::Bool=false) where {VPF,VVPF,MT,VT,DIII,DII,BOO,FLOA,P<:Integer,GRID}
    
    # Extract floating point type from the mesh vertices (Point2{G} where G is the float type)
    FloatType = eltype(eltype(rtm.coarse_mesh).parameters[1]) # Gets G from PolyFace2D{G,T}
    
    # Set defaults based on the mesh's floating point precision
    actual_nudge = nudge === nothing ? 100 * eps(FloatType) : nudge
    actual_rtol = rtol === nothing ? 10 * eps(FloatType) : rtol
    
    if method == :exchange
        exchange_ray_tracing!(rtm, rays_tot, actual_rtol, uncertain, actual_nudge)
    elseif method == :direct
        direct_ray_tracing!(rtm, rays_tot, actual_nudge)
    else
        error("Unknown ray tracing method: $method, must be :exchange or :direct")
    end
end