# Minimal workspace with 2 matrices
mutable struct TwoMatrixWorkspace{P}
    matrix1::Matrix{P}    # B → M (system matrix)
    matrix2::Matrix{P}    # K → S_infty (preserved)
    
    # Working vectors
    work_vec1::Vector{P}
    work_vec2::Vector{P}
    work_vec3::Vector{P}
    work_vec4::Vector{P}
    work_vec5::Vector{P}
    
    # Physical property vectors
    Area::Vector{P}
    epsw::Vector{P}
    Tw::Vector{P}
    qw::Vector{P}
    bin_Qw_known::Vector{Int}
    Volume::Vector{P}
    kappa_g::Vector{P}      # NEW: absorption coefficient per volume
    omega_g::Vector{P}      # NEW: scattering albedo per volume
    Tg::Vector{P}
    qg::Vector{P}
    bin_Qg_known::Vector{Int}
    
    # Constructor
    function TwoMatrixWorkspace{P}(N_surfs::Int, N_vols::Int) where P
        n = N_surfs + N_vols
        new{P}(
            zeros(P, n, n), zeros(P, n, n),  # 2 matrices
            zeros(P, n), zeros(P, n), zeros(P, n), zeros(P, n), zeros(P, n),  # 5 vectors
            zeros(P, N_surfs), zeros(P, N_surfs), zeros(P, N_surfs), zeros(P, N_surfs), zeros(Int, N_surfs),
            zeros(P, N_vols), zeros(P, N_vols), zeros(P, N_vols), zeros(P, N_vols), zeros(P, N_vols), zeros(Int, N_vols)
        )
    end
end

# Simplified workspace for 3D surfaces (no volumes!)
mutable struct SurfaceOnlyWorkspace{P}
    matrix1::Matrix{P}    # B → M (system matrix)
    matrix2::Matrix{P}    # K → S_infty (preserved)
    
    # Working vectors
    work_vec1::Vector{P}
    work_vec2::Vector{P}
    work_vec3::Vector{P}
    work_vec4::Vector{P}
    work_vec5::Vector{P}
    
    # Physical property vectors (surfaces only!)
    Area::Vector{P}
    epsw::Vector{P}
    Tw::Vector{P}
    qw::Vector{P}
    bin_Qw_known::Vector{Int}
    
    # Constructor
    function SurfaceOnlyWorkspace{P}(N_surfs::Int) where P
        new{P}(
            zeros(P, N_surfs, N_surfs), zeros(P, N_surfs, N_surfs),  # 2 matrices
            zeros(P, N_surfs), zeros(P, N_surfs), zeros(P, N_surfs), 
            zeros(P, N_surfs), zeros(P, N_surfs),  # 5 vectors
            zeros(P, N_surfs), zeros(P, N_surfs), zeros(P, N_surfs), 
            zeros(P, N_surfs), zeros(Int, N_surfs)
        )
    end
end

# Populate workspace arrays from mesh - NOW with kappa_g and omega_g precomputation
function populateWorkspace!(ws::TwoMatrixWorkspace{P}, mesh::RayTracingDomain2D, spectral_bin::Int) where P
    surf_count = 0
    vol_count = 0
    
    for i in 1:length(mesh.coarse_mesh)
        for (k, fine_face) in enumerate(mesh.fine_mesh[i])
            vol_count += 1
            ws.Volume[vol_count] = fine_face.volume
            
            # NEW: Store kappa and precompute omega
            ws.kappa_g[vol_count] = fine_face.kappa_g[spectral_bin]
            local_kappa = fine_face.kappa_g[spectral_bin]
            local_sigma_s = fine_face.sigma_s_g[spectral_bin]
            local_beta = local_kappa + local_sigma_s
            ws.omega_g[vol_count] = local_beta > 0.0 ? local_sigma_s / local_beta : 0.0
            
            ws.Tg[vol_count] = fine_face.T_in_g
            ws.qg[vol_count] = fine_face.q_in_g
            ws.bin_Qg_known[vol_count] = fine_face.T_in_g < 0.0 ? 1 : 0

            for (wall_idx, is_solid) in enumerate(fine_face.solidWalls)
                if is_solid
                    surf_count += 1
                    ws.Area[surf_count] = fine_face.area[wall_idx]
                    ws.epsw[surf_count] = fine_face.epsilon[wall_idx][spectral_bin]
                    ws.Tw[surf_count] = fine_face.T_in_w[wall_idx]
                    ws.qw[surf_count] = fine_face.q_in_w[wall_idx]
                    ws.bin_Qw_known[surf_count] = fine_face.T_in_w[wall_idx] < 0.0 ? 1 : 0
                end
            end
        end
    end
end

# Populate workspace from 3D domain (surfaces only!)
function populateWorkspace!(ws::SurfaceOnlyWorkspace{P}, domain::ViewFactorDomain3D, 
                                spectral_bin::Int) where P
    surf_count = 0
    
    for superface in domain.facesMesh
        for subface in superface.subFaces
            surf_count += 1
            ws.Area[surf_count] = subface.area
            
            # Handle spectral or grey epsilon
            if isa(subface.epsilon, Vector)
                ws.epsw[surf_count] = subface.epsilon[spectral_bin]
            else
                ws.epsw[surf_count] = subface.epsilon
            end
            
            ws.Tw[surf_count] = subface.T_in_w
            ws.qw[surf_count] = subface.q_in_w
            ws.bin_Qw_known[surf_count] = subface.T_in_w < 0.0 ? 1 : 0
        end
    end
end

@inline function getValueB(ws::TwoMatrixWorkspace{P}, idx::Int, N_surfs::Int) where {P}
    if idx <= N_surfs
        return 1.0 - ws.epsw[idx]
    else
        return ws.omega_g[idx - N_surfs]
    end
end

# Helper to get B value for surfaces (always just reflectivity)
@inline function getValueB(ws::SurfaceOnlyWorkspace{P}, idx::Int) where {P}
    return 1.0 - ws.epsw[idx]
end