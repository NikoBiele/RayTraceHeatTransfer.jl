function (rtm::RayTracingMeshOptim{VPF,VVPF,MT,VT,DIII,DII,BOO,FLOA,GRID})(gas::GasProperties, rays_tot::P; method::Symbol=:exchange, nudge=Float64(eps(Float32)), rtol=Float64(eps(Float32)), uncertain::Bool=false) where {VPF,VVPF,MT,VT,DIII,DII,BOO,FLOA,P<:Integer,GRID}
    if method == :exchange
        exchange_ray_tracing!(rtm, gas, rays_tot, rtol, uncertain, nudge)
    elseif method == :direct
        direct_ray_tracing!(rtm, gas, rays_tot, nudge)
    else
        error("Unknown ray tracing method: $method, must be :exchange or :direct")
    end
end