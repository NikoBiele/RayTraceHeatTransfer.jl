function steadyStateRigorous(Nx,Ny,N_subs,Area,Volume,FSS,FSG,FGS,FGG,
                                fixWalls,epsw_vec,kappa,maxIter,relTol,
                                Tw_init,Tg_init)
    # this function obtains the steady state temperature distribution
    # in a rigorous (slow) fashion

    sigma = 5.670374419e-08 # [W/m^2-K^4] Stefan-Boltzmann constant

    # create vectors for holding incident energy on all zones
    N_surfs = 2*Nx+2*Ny
    N_vols = Nx*Ny*N_subs
    Gw = zeros(N_surfs)
    Gg = zeros(N_vols)

    # transpose the exchange factor matrices
    FSST = transpose(FSS)
    FSGT = transpose(FSG)
    FGST = transpose(FGS)
    FGGT = transpose(FGG)

    # calculate the A-matrix once
    # A-matrix
    ones_minus_epsw = Diagonal(ones(N_surfs) .- epsw_vec)
    A = I - FSST*ones_minus_epsw

    # Do LU-factorization of the A-matrix
    H = lu!(A)

    # set source terms of walls and volumes
    qw = zeros(N_surfs)
    qg = zeros(N_vols)

    # enter the loop
    iter_count = 0
    Gtot_vec = zeros(maxIter)
    Grelabs = zeros(maxIter)

    # initialize temperatures
    Tw = Tw_init
    Tg = Tg_init
    for i = 1:maxIter
        # calculate emissive powers
        Ew = epsw_vec .* sigma .* Area .* Tw.^4
        Eg = 4 * kappa * sigma * Volume .* Tg.^4

        # calculate the b-vector
        b = [FSST FGST]*[Ew; Eg]

        # solve for the total incident energy on the walls
        # use LU-factorization
        y = H.L\b
        Gw = H.U\y

        # calculate incident energy on volumes
        Gg = [FSGT FGGT]*([Ew; Eg] + [(1 .-epsw_vec).*Gw; zeros(N_vols)])

        # update temperatures of walls if not fixed
        for j = 1:N_surfs
            if fixWalls[j] == true
                Tw[j] = Tw_init[j]
            else
                Tw[j] = ((epsw_vec[j]*Gw[j]+qw[j]*Area[j])/(Area[j]*epsw_vec[j]*sigma))^(1/4)
            end
        end

        # temperature of gas volumes always allowed to change
        Tg = ((Gg + qg .*Volume)./(4*kappa*Volume*sigma)).^(1/4)

        Gtot_vec[i] = sum(Gw) + sum(Gg) # for convergence
        if i > 1
            Grelabs[i] = abs((Gtot_vec[i]-Gtot_vec[i-1])/Gtot_vec[i-1])
            if Grelabs[i] < relTol
                iter_count = i
                break
            end
        end

    end

    return Tw, Tg, iter_count, Grelabs
        
end