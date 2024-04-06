"""
    steadyState(mesh::RayTracingMesh,FSS::Matrix{Float64},FSG::Matrix{Float64},
                        FGS::Matrix{Float64},FGG::Matrix{Float64},
                        epsw_in::Matrix{Float64},gas::GasProperties,
                        maxIter::Int64,relTol::Float64,Tw_in::Matrix{Float64},
                        Tg_in::Vector{Float64},qw_in::Matrix{Float64},qg_in::Vector{Float64})

This function obtains the steady state temperature distribution in a rigorous way.
"""
function steadyState(mesh::RayTracingMesh,FSS::Matrix{Float64},FSG::Matrix{Float64},
                                FGS::Matrix{Float64},FGG::Matrix{Float64},
                                epsw_in::Matrix{Float64},gas::GasProperties,
                                Tw_in::Matrix{Float64},Tg_in::Vector{Float64},
                                qw_in::Matrix{Float64},qg_in::Vector{Float64})

    # first check if the length of input vectors are equal (they should be)
    # give an error if this is true
    if size(Tw_in) != size(qw_in)
        error("The size of Tw and qw are not equal.")
    end
    if size(Tg_in) != size(qg_in)
        error("The size of Tg and qg are not equal.")
    end

    # gas volumes: loop to unpack the volumes, temperatures and sources
    V_count = 0
    Volume = zeros(mesh.N_vols)
    Tg = zeros(mesh.N_vols)
    qg = zeros(mesh.N_vols)
    bin_Qg_known = zeros(mesh.N_vols)
    for k = 1:mesh.N_subs
        # if both a temperature and a source term was specified give warning
        if Tg_in[k] >= 0.0 && qg_in[k] != 0.0
            println("Both a temperature and a source term was specified for the gas in sub-enclosure $k, only the temperature will be used.
                    To specify the source term, set the temperature negative.")
        end
        for m = 1:mesh.Nx
            for n = 1:mesh.Ny
                V_count += 1
                Volume[V_count] = mesh.Volume[k,m,n]
                Tg[V_count] = Tg_in[k]
                qg[V_count] = qg_in[k]
                if Tg_in[k] < 0.0
                    bin_Qg_known[V_count] = 1
                end
            end
        end
    end

    Area = zeros(sum(sum(mesh.solidWalls))*mesh.Nx) # put the right areas into a vector
    wallCount = 0
    for i = 1:mesh.N_subs # for each sub enclosure
        for j = 1:4 # for each wall
            if mesh.solidWalls[i][j] # if wall is solid
                # now found out if everything was specified correctly
                wallCount += 1

                # if both a temperature and a source term was specified give a warning
                if Tw_in[i,j] >= 0.0 && qw_in[i,j] != 0.0
                    println("Both a temperature and a source term was specified for sub-enclosure $i, wall $j, only the temperature will be used.
                            To specify the source term, set the temperature negative.")
                end

                # unpack the area of each wall
                Area[mesh.Nx*(wallCount-1)+1:mesh.Nx*wallCount] = mesh.Area[i,(j-1)*mesh.Nx+1:j*mesh.Nx]
            end
        end
    end
    
    # set the emissivities, temperatures and sources from the input
    wallCount = 0
    epsw = zeros(mesh.N_surfs)
    Tw = zeros(mesh.N_surfs)
    qw = zeros(mesh.N_surfs)
    bin_Qw_known = zeros(mesh.N_surfs)
    for k = 1:mesh.N_subs
        for j = 1:4 # for each wall
            if mesh.solidWalls[k][j] # if wall is solid
                for i = 1:mesh.Nx
                    wallCount += 1
                    epsw[wallCount] = epsw_in[k,j]
                    Tw[wallCount] = Tw_in[k,j]
                    qw[wallCount] = qw_in[k,j]
                    if Tw_in[k,j] < 0.0
                        bin_Qw_known[wallCount] = 1
                    end
                end
            end
        end
    end

    sigma = 5.670374419e-08 # [W/m^2-K^4] Stefan-Boltzmann constant

    # transpose the exchange factor matrices
    FSST = transpose(FSS)
    FSGT = transpose(FSG)
    FGST = transpose(FGS)
    FGGT = transpose(FGG)
    FT = [FSST FGST; FSGT FGGT] # the full transposed exchange factor matrix

    # construct some matrices used for solving
    I_Q_known = Diagonal([bin_Qw_known; bin_Qg_known])
    ones_minus_epsw = ones(mesh.N_surfs) .- epsw
    zeros_gas = zeros(mesh.N_vols)
    I_ones_minus_epsw = Diagonal([ones_minus_epsw; zeros_gas])
    I_epsw_and_ones = Diagonal([epsw; ones(mesh.N_vols)])

    # found out which emissive powers and sources are known
    # walls first
    Ew_known = zeros(mesh.N_surfs)
    Qw_known = zeros(mesh.N_surfs)
    for i = 1:mesh.N_surfs
        if Tw[i] >= 0.0
            # if temperature is known, calculate emissive power
            Ew_known[i] = epsw[i]*sigma*Area[i]*Tw[i]^4
        else
            # if temperature is unknown, calculate source
            Qw_known[i] = qw[i]*Area[i]
        end
    end
    # volumes next
    Eg_known = zeros(mesh.N_vols)
    Qg_known = zeros(mesh.N_vols)
    for i = 1:mesh.N_vols
        if Tg[i] >= 0.0
            # if temperature is known, calculate emissive power
            Eg_known[i] = 4*gas.kappa*sigma*Volume[i]*Tg[i]^4
        else
            # if temperature is unknown, calculate source
            Qg_known[i] = qg[i]*Volume[i]
        end
    end
    E_known = [Ew_known; Eg_known] # all known emissive powers
    Q_known = [Qw_known; Qg_known] # all known sources

    # setup and solve the rigorous large system
    A_matrix = I - FT * (I_epsw_and_ones*I_Q_known + I_ones_minus_epsw)
    b_vector = FT * (E_known .+ Q_known)
    G = A_matrix\b_vector

    # calculate unknown emissive powers (temperatures) and sources
    # walls first
    for i = 1:mesh.N_surfs
        Ew = epsw[i]*G[i] + Q_known[i]
        if Tw[i] < 0.0
            # if temperature is unknown, calculate it
            Tw[i] = (Ew/(epsw[i]*sigma*Area[i]))^(1/4)
        else
            # if temperature is known, calculate source
            Qw = Ew - epsw[i]*G[i]
            qw[i] = Qw/Area[i]
        end
    end
    # volumes next
    for i = mesh.N_surfs+1:mesh.N_surfs+mesh.N_vols
        Eg = G[i] + Q_known[i]
        if Tg[i-mesh.N_surfs] < 0.0
            # if temperature is unknown, calculate it
            Tg[i-mesh.N_surfs] = (Eg/(4*gas.kappa*sigma*Volume[i-mesh.N_surfs]))^(1/4)
        else
            # if temperature is known, calculate source
            Qg = Eg - G[i]
            qg[i-mesh.N_surfs] = Qg/Volume[i-mesh.N_surfs]
        end
    end

    return Tw, Tg, qw, qg
        
end