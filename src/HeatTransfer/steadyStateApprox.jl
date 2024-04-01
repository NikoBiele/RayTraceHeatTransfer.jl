"""
    steadyStateRigorous(mesh::RayTracingMesh,FSS::Matrix{Float64},FSG::Matrix{Float64},
                        FGS::Matrix{Float64},FGG::Matrix{Float64},
                        epsw_in::Matrix{Float64},gas::GasProperties,
                        maxIter::Int64,relTol::Float64,Tw_in::Matrix{Float64},
                        Tg_in::Vector{Float64},qw_in::Matrix{Float64},qg_in::Vector{Float64})

This function obtains the steady state temperature distribution in a rigorous way.
"""
function steadyStateApprox(mesh::RayTracingMesh,FSS::Matrix{Float64},FSG::Matrix{Float64},
                                FGS::Matrix{Float64},FGG::Matrix{Float64},
                                epsw_in::Matrix{Float64},gas::GasProperties,
                                maxIter::Int64,relTol::Float64,Tw_in::Matrix{Float64},
                                Tg_in::Vector{Float64},qw_in::Matrix{Float64},qg_in::Vector{Float64})

    # first check if the length of input vectors are equal (they should be)
    # give an error if this is true
    if size(Tw_in) != size(qw_in)
        error("The size of Tw and qw are not equal.")
    end
    if size(Tg_in) != size(qg_in)
        error("The size of Tg and qg are not equal.")
    end

    # found out which wall temperatures and source terms are fixed by the user
    fixWallTemp = Vector{Bool}(undef,mesh.N_surfs)
    
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

                # if the wall temperature was specified the wall temperature is fixed
                # then the corresponding source term is set to zero
                if Tw_in[i,j] < 0
                    fixWallTemp[mesh.Nx*(wallCount-1)+1:mesh.Nx*wallCount] .= false
                else
                    fixWallTemp[mesh.Nx*(wallCount-1)+1:mesh.Nx*wallCount] .= true
                end

                # unpack the area of each wall
                Area[mesh.Nx*(wallCount-1)+1:mesh.Nx*wallCount] = mesh.Area[i,(j-1)*mesh.Nx+1:j*mesh.Nx]
            end
        end
    end

    # found out which gas temperatures and source terms are fixed by the user
    fixGasTemp = Vector{Bool}(undef,mesh.N_subs*mesh.Nx*mesh.Ny)
    for k = 1:mesh.N_subs # for each sub enclosure
        # if both a temperature and a source term was specified give warning
        if Tg_in[k] >= 0.0 && qg_in[k] != 0.0
            println("Both a temperature and a source term was specified for the gas in sub-enclosure $k, only the temperature will be used.
                    To specify the source term, set the temperature negative.")
        end
        
        # if the wall temperature was specified the wall temperature is fixed
        # then the corresponding source term is set to zero
        if Tg_in[k] < 0
            fixGasTemp[mesh.Nx*mesh.Ny*(k-1)+1:mesh.Nx*mesh.Ny*k] .= false
        else
            fixGasTemp[mesh.Nx*mesh.Ny*(k-1)+1:mesh.Nx*mesh.Ny*k] .= true
        end
    end

    # loop to unpack the volumes, temperatures and sources
    V_count = 0
    Volume = zeros(mesh.N_subs*mesh.Nx*mesh.Ny)
    Tg = zeros(mesh.N_subs*mesh.Nx*mesh.Ny)
    qg = zeros(mesh.N_subs*mesh.Nx*mesh.Ny)
    for k = 1:mesh.N_subs
        for m = 1:mesh.Nx
            for n = 1:mesh.Ny
                V_count += 1
                Volume[V_count] = mesh.Volume[k,m,n]
                Tg[V_count] = Tg_in[k]
                qg[V_count] = qg_in[k]
            end
        end
    end

    # set the emissivities, temperatures and sources from the input
    wallCount = 0
    epsw = zeros(mesh.N_surfs)
    Tw = zeros(mesh.N_surfs)
    qw = zeros(mesh.N_surfs)
    for k = 1:mesh.N_subs
        for j = 1:4 # for each wall
            if mesh.solidWalls[k][j] # if wall is solid
                for i = 1:mesh.Nx
                    wallCount += 1
                    epsw[wallCount] = epsw_in[k,j]
                    Tw[wallCount] = Tw_in[k,j]
                    qw[wallCount] = qw_in[k,j]
                end
            end
        end
    end

    sigma = 5.670374419e-08 # [W/m^2-K^4] Stefan-Boltzmann constant

    # create vectors for holding incident energy on all zones
    Gw = zeros(mesh.N_surfs)
    Gg = zeros(mesh.N_vols)

    # transpose the exchange factor matrices
    FSST = transpose(FSS)
    FSGT = transpose(FSG)
    FGST = transpose(FGS)
    FGGT = transpose(FGG)

    # calculate the A-matrix once
    ones_minus_epsw = Diagonal(ones(mesh.N_surfs) .- epsw)
    A = I - FSST*ones_minus_epsw

    # Do LU-factorization of the A-matrix
    H = lu!(A)

    # enter the loop
    iter_count = 0
    Gtot_vec = zeros(maxIter)
    Grelabs = zeros(maxIter)
    T_save = zeros(length(Tw)+length(Tg),maxIter)
    T_save[1:length(Tw),1] = Tw
    T_save[length(Tw)+1:length(Tw)+length(Tg),1] = Tg

    for i = 2:maxIter
        # calculate emissive powers
        Ew = epsw .* sigma .* Area .* Tw.^4
        Eg = 4 * gas.kappa * sigma * Volume .* Tg.^4

        # calculate the b-vector
        b = [FSST FGST]*[Ew; Eg]

        # solve for the total incident energy on the walls
        # use LU-factorization
        y = H.L\b
        Gw = H.U\y

        # calculate incident energy on volumes
        Gg = [FSGT FGGT]*([Ew; Eg] + [(1 .-epsw).*Gw; zeros(mesh.N_vols)])

        # update temperatures of walls if not fixed
        for j = 1:mesh.N_surfs
            if fixWallTemp[j] == false
                # only update wall temperatures which has not been specified
                Tw[j] = ((epsw[j]*Gw[j]+qw[j]*Area[j])/(Area[j]*epsw[j]*sigma))^(1/4)
            end
        end

        # temperature of gas volumes always allowed to change
        for j = 1:mesh.N_vols
            if fixGasTemp[j] == false
                Tg[j] = ((Gg[j] + qg[j]*Volume[j])./(4*gas.kappa*Volume[j]*sigma)).^(1/4)
            end
        end

        Gtot_vec[i] = sum(Gw) + sum(Gg) # for convergence
        if i > 1
            Grelabs[i] = abs((Gtot_vec[i]-Gtot_vec[i-1])/Gtot_vec[i-1])
            if Grelabs[i] < relTol
                iter_count = i
                break
            end
        end
    end

    # return temperatures, fluxes and iteration information
    return Tw, Tg, Gw, Gg, iter_count, Grelabs
        
end