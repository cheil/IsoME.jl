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


"""
    calc_Ne_Sc(mu, Ne_nsc, itemp, wsi, dos_en, dos, znormip, deltaip, shiftip)

Calculate the number of electrons in the sc state for a given chemical
potential minus the number of electrons in the normal state
According to Lucrezi, Communication Physics, (2024) 7:33, eq. (16)
or Lee, Computational Materials (2023) 9:156, eq. (32) 
(See also Overleaf/Matsubara_sums)
"""
function calc_Ne_Sc(mu, Ne_nsc, itemp, wsi, dos_en, dos, znormip, deltaip, shiftip)
    # eq (9) & (11) overleaf
    diff = dos_en .- mu
    theta = (wsi' .* znormip') .^ 2 .+ (diff .+ shiftip') .^ 2 .+ (znormip' .* deltaip) .^ 2
    summand_sc = (diff .+ shiftip') ./ theta .- diff ./ (wsi'.^ 2 .+ diff.^ 2)
    
    # matsubara sum
    summand_sc = dropdims((sum(summand_sc, dims=2)), dims=2)

    summand_sc = 2 * fermiFcn(dos_en, mu, itemp) - 4 * kb * itemp * summand_sc

    # diff between Ne in normal and sc state
    Ne_sc =  Ne_nsc - trapz(dos_en, dos.*summand_sc) 

    return Ne_sc

end


"""
    update_mu_own(itemp, wsi, ef, dos_en, dos, znormip, deltaip, shiftip)

Routine to update chemical potential s.t. the number of electrons stays fixed

The routine uses the bisection method to find a value for the chemical potential
where the amount of electrons in the SC state is equal to the normal state
"""
function update_mu_own(itemp, wsi, dos_en, dos, znormip, deltaip, shiftip, idxEncut, outdir)

    # delta as row vector, needed if no weep
    if size(deltaip, 2) == 1
        deltaip = deltaip'
    else
        deltaip = deltaip[idxEncut[1]:idxEncut[2],:]
    end

    ### Calculate N_e in the non-SC state
    Ne_nsc = trapz(dos_en[idxEncut[1]:idxEncut[2]], 2 .* fermiFcn(dos_en[idxEncut[1]:idxEncut[2]], 0.0, itemp) .* dos[idxEncut[1]:idxEncut[2]])   

    # call calc_Ne_Sc with first argument unspecified
    fmu(x) = calc_Ne_Sc(x, Ne_nsc, itemp, wsi, dos_en[idxEncut[1]:idxEncut[2]], dos[idxEncut[1]:idxEncut[2]], znormip, deltaip, shiftip)  

    ### starting values for mu
    mu0 = -100
    mu1 = +100
    fmu0 = fmu(mu0)
    fmu1 = fmu(mu1)

    ### wrong slope 
    if fmu1 > fmu0
        # plot electron number vs mu
        mu_error = range(mu0 - 100, mu1 + 100, 200)
        Ne_error = zeros(size(mu_error))
        for k in eachindex(mu_error)
            Ne_error[k] = fmu(mu_error[k])
        end

        p = plot(mu_error, Ne_error, label="Ne_nsc - Ne_sc", title="Ne in normal state minus sc state")
        #vline(p, [mu0, mu1], label="mu")
        savefig(outdir*"muError.png")

        error("The number of electrons decreases with increasing mu! Try a larger omega_c")
    end

    ### find minimum interval around ef in which a sign change occurs
    iter = 0
    mu0error = mu0
    mu1error = mu1
    while fmu0 * fmu1 > 0
        if sign(fmu0) < 0
            mu0error -= 100
            mu0 -= 100
            mu1 -= 100
            fmu0 = fmu(mu0)
        else
            mu0 += 100
            mu1 += 100
            mu1error += 100
            fmu1 = fmu(mu1)
        end
        iter += 1
        if iter > 200
            error("Error in mu update - Couldn't find a root. Please check your input files, in particular the dos-file.")
        end
    end

    ### calc new mu using the bisection method
    mu = bisection(fmu, mu0, mu1)
    #mu = find_zero(fmu, [mu0, mu1])
    

    return mu


end

