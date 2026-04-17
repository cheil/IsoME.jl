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
function diff_Ne(mu, Ne_nsc, itemp, wsi, dos_en, dos, znormip, deltaip, shiftip)
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
    root_finding(fmu)

Find the root of f(mu) = Ne_sc(mu) - Ne_nsc = 0.
"""
function root_finding(fmu, outdir)

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

        error("The number of electrons decreases with increasing mu!")
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


"""
    update_mu_own(itemp, wsi, ef, dos_en, dos, znormip, deltaip, shiftip)

Routine to update chemical potential s.t. the number of electrons stays fixed

The routine uses the bisection method to find a value for the chemical potential
where the amount of electrons in the SC state is equal to the normal state
"""
function update_mu_own(itemp, wsi, dos_en, dos, znormip, deltaip, shiftip, idxShiftcut, outdir)

    # delta as row vector, needed if no weep
    if size(deltaip, 2) == 1
        deltaip = deltaip'
    else
        deltaip = deltaip[idxShiftcut[1]:idxShiftcut[2],:]
    end

    ### Calculate N_e in the non-SC state
    Ne_nsc = trapz(dos_en[idxShiftcut[1]:idxShiftcut[2]], 2 .* fermiFcn(dos_en[idxShiftcut[1]:idxShiftcut[2]], 0.0, itemp) .* dos[idxShiftcut[1]:idxShiftcut[2]])   

    # call calc_Ne_Sc with first argument unspecified
    fmu(x) = diff_Ne(x, Ne_nsc, itemp, wsi, dos_en[idxShiftcut[1]:idxShiftcut[2]], dos[idxShiftcut[1]:idxShiftcut[2]], znormip, deltaip, shiftip)  

    mu = root_finding(fmu, outdir)
    
    return mu

end



#########################################################
# --------------------- Real Axis --------------------- #
#########################################################
"""
    mu_update_real_axis()

Update the chemical potential to conserve charge neutrality - real axis implementation.
"""
function mu_update_real_axis(itemp, w_static, w_static_chi, w_prime, dos_en, dos, znormip, deltaip, shiftip, outdir, FLIP, i_it)

    # idx_shiftcut in dos ???

    ### Calculate N_e in the non-SC state
    Ne_nsc = 2 .* trapz(dos_en, fermiFcn(dos_en, 0.0, itemp) .* dos)   

    # call calc_Ne_Sc with first argument unspecified
    fmu(x) = diff_Ne_realAxis(x, Ne_nsc, itemp, w_static, w_static_chi, w_prime, dos_en, dos, znormip, deltaip, shiftip, FLIP, outdir, i_it)

    fmu(-4000)

    ytest = Vector{Float64}()
    for xtest in range(-5000, 0,100)
        push!(ytest, fmu(xtest))
    end
    #println(ytest)
    plot(range(-5000, 0,100), ytest)
    savefig("mutest.png")

    mu = find_zero(fmu, 0.0, Order1())

    #mu = root_finding(fmu, outdir)
    
    return mu
end


"""
    diff_Ne_realAxis(mu, Ne_nsc, itemp, w_prime, dos_en, dos, znormip, deltaip, shiftip)

Calculate the number of electrons in the sc state for a given chemical
potential minus the number of electrons in the normal state
"""
function diff_Ne_realAxis(mu, Ne_nsc, itemp, w_static, w_static_chi, w_prime, dos_en, dos, znormip, deltaip, shiftip, FLIP, outdir, i_it)

    phiphip = deltaip.*znormip

    # interpolate Z,χ,ϕ onto ω'-integration grid
    Z_itp = linear_interpolation(w_static, znormip, extrapolation_bc=Flat())
    phi_itp = linear_interpolation(w_static, phiphip, extrapolation_bc=Flat())
    shift_itp = linear_interpolation(w_static_chi, shiftip, extrapolation_bc=Flat())

    Z_ongrid = Z_itp.(w_prime)
    phi_ongrid = phi_itp.(w_prime)
    shift_ongrid = shift_itp.(w_prime)

    # --------------- Causality --------------- #
    M0, M1, ε_p, Rplus, Rminus, Iplus, Iminus, I0_pp, I1_pp, I2_pp, I3_pp, I0_mm, I1_mm, I2_mm, I3_mm = epsilon_helpers(dos_en, dos, Z_ongrid, phi_ongrid, shift_ongrid .- mu, w_prime)


    #println("I0_pp", I0_pp[1:10])

    z_integrand = eval_spectral_integrals(w_prime.*Z_ongrid, M0, M1, ε_p, Rplus, Rminus, Iplus, Iminus, I0_pp, I0_mm, I1_pp, I1_mm, I2_pp, I2_mm, I3_pp, I3_mm, false, python_style=true, epsilon=dos_en)     # ε + (χ(ω) - μ_F)
    FLIP = ifelse.(-abs.(z_integrand) .== z_integrand, 1, -1) 

    #print(z_integrand[1:10])

    #error("s")


    # ------------- ε-integration ------------- #
    omega_integrand = eval_spectral_integrals(shift_ongrid-mu, M0, M1, ε_p, Rplus, Rminus, Iplus, Iminus, I0_pp, I0_mm, I1_pp, I1_mm, I2_pp, I2_mm, I3_pp, I3_mm, true)     # ε + (χ(ω) - μ_F)
    
    # plot(w_prime, omega_integrand)
    # savefig("integrand_"*string(i_it)*".png")

    # error("a")

    omega_integrand .*= FLIP
    β = 1/(kb*itemp)
    omega_integrand .*= tanh.(β*w_prime/2)

    # ------------- ω-integration ------------- #
    Ne_sc = trapz(dos_en, dos) + 2/π*trapz(w_prime, omega_integrand)

    # diff between Ne in normal and sc state
    root_eq =  Ne_nsc - Ne_sc


    return root_eq
end