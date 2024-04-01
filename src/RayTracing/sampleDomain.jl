"""
    sampleDomain(mesh::RayTracingMesh,gas::GasProperties,
                N_rays::Int64,nthreads::Int64,displayWhileTracing::Bool)

This function ray traces the entire domain, both surfaces and volumes.
"""
function sampleDomain(mesh::RayTracingMesh,gas::GasProperties,
                    N_rays::Int64,nthreads::Int64,displayWhileTracing::Bool)

    ### SAMPLE SURFACES
    FSS, FSG = sampleSurfaces(mesh,gas,N_rays,nthreads,displayWhileTracing);

    ### SAMPLE VOLUMES
    FGS, FGG = sampleVolumes(mesh,gas,N_rays,nthreads,displayWhileTracing);

    return FSS, FSG, FGS, FGG
end