"""
    File containing an update routine to find the chemical potential 
    which keeps the number of electrons constant in the sc phase
    
    Julia Packages:
        - 

    Comments:
        - 

"""



function fermiFcn(epsilon, mu, T)
    """
    Fermi-Dirac distribution

    --------------------------------------------------------------------
    Input:
        mu:         chemical potential/fermi energy
        epsilon:    energies 
        T:          Temperature

    --------------------------------------------------------------------
    Output:
        nF:         Fermi-Dirac distribution

    --------------------------------------------------------------------
    Comments:
        -
        
    -------------------------------------------------------------------- 
    """


    nF = 1.0 ./(exp.((epsilon.-mu)/(kb*T)).+1)

    return nF

end


function calc_Ne_Sc(mu, Ne_nsc, itemp, wsi, M, dos_en, dos, znormip, deltaip, shiftip)
    """
    Calculate the number of electrons in the sc state for a given chemical
    potential minus the number of electrons in the normal state
    According to Lucrezi, Communication Physics, (2024) 7:33, eq. (16)
    or Lee, Computational Materials (2023) 9:156, eq. (32) 
    (See also Overleaf/Matsubara_sums)

    --------------------------------------------------------------------
    Input:
        mu:         chemical potential
        Ne_nsc:     Number of electrons in normal state
        itemp:      Temperature 
        wsi:        vector containing matsubara frequencies
        M:          number of highest matsubara frequency
        dos_en:     energy grid points
        dos:        density of states
        znormip:    Z of previous iteration
        deltaip:    Gap of previous iteration
        shiftip:    shift of previous iteration

    --------------------------------------------------------------------
    Output:
        Ne:         Number of electrons in sc state

    --------------------------------------------------------------------
    Comments:
        -
        
    -------------------------------------------------------------------- 
    """

    # eq (9) & (11) overleaf
    delta = dos_en .- mu
    theta = (wsi' .* znormip') .^ 2 .+ (delta .+ shiftip') .^ 2 .+ (znormip' .* deltaip) .^ 2
    integrand_Nsc = (delta .+ shiftip') ./ theta .- delta ./ (wsi' .^ 2 .+ delta .^ 2)
    
    # matsubara sum
    integrand_Nsc = dropdims((sum(integrand_Nsc, dims=2)), dims=2)
    integrand_Nsc = 2 * fermiFcn(dos_en, mu, itemp) - 4 * kb * itemp * integrand_Nsc

    # diff between Ne in normal and sc state
    Ne_sc =  Ne_nsc - trapz(dos_en, dos.*integrand_Nsc) 

    return Ne_sc

end



function update_mu_own(itemp, wsi, M, muintr, dos_en, dos, znormip, deltaip, shiftip)
    """
    Routine to update chemical potential (own version not relying on nel)
    Newton - Raphson, doesn't work here, it diverges after a few iterations
    I'm using the secant method that seems to work very well

    --------------------------------------------------------------------
    Input:
        itemp:      Temperature 
        wsi:        vector containing matsubara frequencies
        M:          number of highest matsubara frequency
        muintr:
        ndos:       number of energy grid points
        dos_en:     energy grid points
        dos:        density of states
        znormip:    Z of previous iteration
        deltaip:    Gap of previous iteration
        shiftip:    shift of previous iteration

    --------------------------------------------------------------------
    Output:
        mu:         updated chemical potential

    --------------------------------------------------------------------
    Comments:
        -
        
    -------------------------------------------------------------------- 
    """



    ### Intialize
    # max iterations
    Nit = 5000

    # length energies
    ndos = size(dos_en, 1)

    ### Calculate f(mu) in the non-SC state
    integrand_fmu_nsc = zeros(ndos)
    for iw in 1:1:M+1
        delta = dos_en .- ef
        theta = (wsi[iw]) .^ 2 .+ (dos_en .- ef) .^ 2     

        integrand_fmu_nsc = integrand_fmu_nsc .+ delta ./ theta
    end
    Ne_nsc = trapz(dos_en, 2 .* fermiFcn(dos_en, ef, itemp) .* dos)

    # delta as row vector, needed if no weep
    if include_Weep == 0
        deltaip = deltaip'
    end

    # call calc_Ne_Sc with first argument unspecified
    fmu(x) = calc_Ne_Sc(x, Ne_nsc, itemp, wsi, M, dos_en, dos, znormip, deltaip, shiftip)  

    ### starting values for mu
    mu0 = ef - 100
    mu1 = ef + 100
    fmu0 = fmu(mu0)
    fmu1 = fmu(mu1)

    ### wrong slope 
    if fmu1 > fmu0
        # plot electron number vs mu
        mu_error = range(mu0 - 500, mu1 + 400, 100)
        Ne_error = zeros(size(mu_error))
        for k in eachindex(mu_error)
            Ne_error[k] = fmu(mu_error[k])
        end

        p = plot(mu_error, Ne_error, label="Ne_nsc - Ne_sc", title="Ne in normal state minus sc state")
        #vline(p, [mu0, mu1], label="mu")
        savefig("muError.png")

        error("The number of electrons decreases with increasing mu! Try a larger omega_c")
    end

    ### find minimum interval around ef in which a sign change occurs
    iter = 0
    while fmu0 * fmu1 > 0
        if sign(fmu0) < 0
            mu0 -= 100
            mu1 -= 100
            fmu0 = fmu(mu0)
        else
            mu0 += 100
            mu1 += 100
            fmu1 = fmu(mu1)
        end
        iter += 1
        if iter > 200
            error("Error in mu update. Couldn't find a root within an interval of " * string(mu0) * " meV and " * string(mu1) * " meV")
        end
    end

    ### calc new mu using the bisection method    
    mu = bisection(fmu, mu0, mu1)
    #mu = RegulaFalsi(fmu, mu0, mu1)
    #mu3 = find_zero(fmu, [mu0, mu1])

    #println(mu)
    return mu


end



