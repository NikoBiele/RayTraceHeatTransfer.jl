## Public and private types and functions

### Public functions and types for geometry creation

```@meta
CurrentModule = RayTraceHeatTransfer
```

```@docs
Geometry.SubEnclosure
```

```@docs
Geometry.displayGeometry(geometry::Vector{SubEnclosure})
```

```@docs
Geometry.meshGeometry(subs::Vector{SubEnclosure},Ny::Int64,Nx::Int64) 
```

```@docs
Geometry.RayTracingMesh
```

```@docs
Geometry.displayMesh(mesh::RayTracingMesh)
```

```@docs
Geometry.GasProperties
```

### Public functions for ray tracing

```@docs
RayTracing.sampleDomain(mesh::RayTracingMesh,gas::GasProperties,
                N_rays::Int64,nthreads::Int64,displayWhileTracing::Bool)
```

```@docs
RayTracing.writeMatricesToCSV(FSS::Matrix{Float64},FSG::Matrix{Float64},
                        FGS::Matrix{Float64},FGG::Matrix{Float64},
                        mesh::RayTracingMesh,gas::GasProperties,
                        N_rays::Int64)
```

### Public functions for heat transfer calculation

```@docs
HeatTransfer.readMatricesFromCSV(N_subs::Int64,Ndim::Int64,kappa::Float64,
                            sigma_s::Float64,N_rays::Int64)
```

```@docs
HeatTransfer.steadyState(mesh::RayTracingMesh,FSS::Matrix{Float64},FSG::Matrix{Float64},
                    FGS::Matrix{Float64},FGG::Matrix{Float64},
                    epsw_in::Matrix{Float64},gas::GasProperties,
                    Tw_in::Matrix{Float64},Tg_in::Vector{Float64},
                    qw_in::Matrix{Float64},qg_in::Vector{Float64})
```

```@docs
HeatTransfer.plotTemperatureField(mesh::RayTracingMesh,Tg::Vector{Float64},Tw=nothing)
```

```@docs
HeatTransfer.validCrosbieSchrenker(N_rays_tot::Int64,Ndim::Int64,Tw_hot::Float64)
```

### Public utility 3D view factor function

```@docs
ViewFactor3D.viewFactor(POLY_A::Matrix{Float64}, POLY_B::Matrix{Float64})
```

### Private functions

#### Private geometry functions

```@docs
Geometry.areaVolumeMesh(Nx::Int64,Ny::Int64,N_subs::Int64,
            point1::Matrix{SVector{2,Float64}},point2::Matrix{SVector{2,Float64}},
            point3::Matrix{SVector{2,Float64}},point4::Matrix{SVector{2,Float64}})
```

```@docs
Geometry.defineSolidWalls(subs::Vector{SubEnclosure})
```

```@docs
Geometry.localWalls(point1::Vector{Float64},point2::Vector{Float64},
                point3::Vector{Float64},point4::Vector{Float64})
```

#### Private functions for ray tracing

```@docs
RayTracing.coefs_exchange(mesh::RayTracingMesh, Wall_absorbX::Array{Float64, 3},
                        Wall_absorbY::Array{Float64, 3},
                        N_abs_gas::Array{Float64, 3},RayCountTotal::Float64)
```

```@docs
RayTracing.distToSurface(point::SVector{2,Float64}, i1::SVector{2,Float64},
                    mesh::RayTracingMesh, N_subs_count::Int64)
```

```@docs
RayTracing.isotropicScatter()
```

```@docs
RayTracing.lambertSample3D()
```

```@docs
RayTracing.rayTracing_CPU(subNumber::Int64,wallNumber::Int64,
                        sampleLeftRight::Bool,sampleTopBottom::Bool,
                        mesh::RayTracingMesh, gas::GasProperties, N_rays::Int64,
                        wallEmission::Bool,volumeEmission::Bool,
                        displayWhileTracing::Bool,nthreads::Int64,
                        xCountSample::Int64, yCountSample::Int64)
```

```@docs
RayTracing.sampleSurface(mesh::RayTracingMesh, subNumber::Int64, wallNumber::Int64,
                            xCountSample::Int64, yCountSample::Int64,
                            sampleLeftRight::Bool, sampleTopBottom::Bool)
```

```@docs
RayTracing.sampleSurfaces(mesh::RayTracingMesh,gas::GasProperties,
                    N_rays::Int64,nthreads::Int64,displayWhileTracing::Bool)
```

```@docs
RayTracing.sampleVolume(mesh::RayTracingMesh, xCount::Int64, yCount::Int64)
```

```@docs
RayTracing.sampleVolumes(mesh::RayTracingMesh,gas::GasProperties,
                    N_rays::Int64,nthreads::Int64,displayWhileTracing::Bool)
```

```@docs
RayTracing.solidWalls(mesh::RayTracingMesh, N_subs_count::Int64)
```

```@docs
RayTracing.whichCell(point::SVector{2,Float64},mesh::RayTracingMesh)
```

```@docs
RayTracing.whichSubEnclosure(point::SVector{2,Float64},mesh::RayTracingMesh)
```

#### Private functions for view factor function

```@docs
ViewFactor3D.Cl(theta::Float64, almostZero::Float64)
```

```@docs
ViewFactor3D.edgePairParameters(Po::Matrix{Float64}, Pf::Matrix{Float64},
                        Qo::Matrix{Float64}, Qf::Matrix{Float64},
                        almostZero::Float64)
```

```@docs
ViewFactor3D.f(s::Float64, l::Float64, alpha::Float64,
        cosAlpha::Float64, sinAlpha::Float64,
        d::Float64, almostZero::Float64)
```

```@docs
ViewFactor3D.fParallel(s::Float64, l::Float64,
            d::Float64, almostZero::Float64)
```

```@docs
ViewFactor3D.imagLi_2(mag::Float64, angle::Float64,
                    almostZero::Float64)
```

## Index

```@index
```