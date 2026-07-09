"""
    File containing an update routine to find the chemical potential 
    which keeps the number of electrons constant in the sc phase
    
    Julia Packages:
        - 

    Comments:
        - 

"""


"""
    fermiFcn(epsilon, mu, T)

Fermi-Dirac distribution n_F(epsilon-mu; T)
"""
function fermiFcn(epsilon, mu, T)

    nF = 1.0 ./(exp.((epsilon.-mu)/(kb*T)).+1)

    return nF

end


"""
    root_finding(fmu)

Find the root of f(mu) = Ne_nsc(mu) - Ne_sc = 0.
"""
function root_finding(fmu, outdir, fermi_level, i_it=0)

    ### starting values for mu
    mu0 = fermi_level-100
    mu1 = fermi_level+100
    fmu0 = fmu(mu0)
    fmu1 = fmu(mu1)

    ### wrong slope 
    if fmu1 > fmu0
        # plot electron number vs mu
        mu_error = range(mu0 - 100, mu1 + 100, 50)
        Ne_error = zeros(size(mu_error))
        for k in eachindex(mu_error)
            Ne_error[k] = fmu(mu_error[k])
        end

        p = plot(mu_error, Ne_error, label="Ne_nsc - Ne_sc", title="Ne in normal state minus sc state")
        #vline(p, [mu0, mu1], label="mu")
        savefig(outdir*"muError.png")
        savePlotData(outdir*"muError.dat", "#  μ / meV       Ne_nsc - Ne_sc", mu_error, Ne_error)

        error("The number of electrons decreases with increasing mu!")
    end

    ### find minimum interval around ef in which a sign change occurs
    iter = 0
    mu_error = [mu0, mu1]
    fmu_error = [fmu0, fmu1]
    while fmu0 * fmu1 > 0
        if sign(fmu0) < 0
            mu0 -= 20
            mu1 -= 20
            fmu1 = fmu0
            fmu0 = fmu(mu0)
            pushfirst!(mu_error, mu0)
            pushfirst!(fmu_error, fmu0)
        else
            mu0 += 20
            mu1 += 20
            fmu0 = fmu1
            fmu1 = fmu(mu1)
            push!(mu_error, mu1)
            push!(fmu_error, fmu1)
        end
        iter += 1
        if iter > 100    # 2 eV
            plot(mu_error, fmu_error, label="Ne_nsc - Ne_sc", title="Ne in normal state minus sc state")
            savefig(outdir*"muError.png")

            mu0error = mu_error[1]
            mu1error = mu_error[end]
            error("Error in mu update - Couldn't find a root in the interval [$mu0error,$mu1error]. Please check your input files, in particular the dos-file.")
        end
    end

    ### calc new mu using the RegulaFalsi method
    mu = RegulaFalsi(fmu, mu0, mu1, 1e-3, 1e-6)

    return mu
end



############################################################
# -------------------- Matsubara axis -------------------- #
############################################################
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
    update_mu_own(itemp, wsi, ef, dos_en, dos, znormip, deltaip, shiftip)

Routine to update chemical potential to fix the number of electrons
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

    mu = root_finding(fmu, outdir, 0)
    
    return mu

end



############################################################
# ---------------------- Real Axis ----------------------- #
############################################################
"""
    mu_update_real_axis()

Update the chemical potential to conserve charge neutrality - real axis implementation.
"""
function mu_update_real_axis(itemp, fermi_level, w_static, w_static_chi, w_prime, dos_en, dos, znormip, phiphip, shiftip, outdir, i_it)

    ### Calculate N_e in the non-SC state
    Ne_nsc = 2 .* trapz(dos_en, fermiFcn(dos_en, 0.0, itemp) .* dos)   

    # interpolate Z,χ,ϕ 
    Z_itp = linear_interpolation(w_static, znormip, extrapolation_bc=Flat())
    phi_itp = linear_interpolation(w_static, phiphip, extrapolation_bc=Flat())
    shift_itp = linear_interpolation(w_static_chi, shiftip, extrapolation_bc=Flat())

    Z_ongrid = Z_itp.(w_prime)
    phi_ongrid = phi_itp.(w_prime)
    shift_ongrid = shift_itp.(w_prime)

    # piecewise-linear DOS segments
    eps_1, eps_2 = transpose(dos_en[1:end-1]), transpose(dos_en[2:end])
    dos_1, dos_2 = transpose(dos[1:end-1]), transpose(dos[2:end])
    electronic_spec = (eps_1, eps_2, dos_1, dos_2, eps_2 .- eps_1, dos_2 .- dos_1)

    dos_int = trapz(dos_en, dos)

    # temperature
    β = 1/(kb*itemp)
    tanhw = tanh.(β .* w_prime ./ 2)
    
    # call calc_Ne_Sc with first argument unspecified
    fmu(x)  = diff_Ne_realAxis(x, Ne_nsc, tanhw, w_prime, electronic_spec, dos_int, Z_ongrid, phi_ongrid, shift_ongrid)

    mu = root_finding(fmu, outdir, fermi_level, i_it)
    #@time mu2 = find_zero((fmu, dfmu), fermi_level, Roots.LithBoonkkampIJzerman(3, 1))
   
    return mu
end


"""
    diff_Ne_realAxis(mu, Ne_nsc, itemp, w_prime, dos_en, dos, znormip, deltaip, shiftip)

Calculate the number of electrons in the sc state for a given chemical
potential minus the number of electrons in the normal state
"""
function diff_Ne_realAxis(mu, Ne_nsc, tanhw, w_prime, electronic_spec, dos_int, Z_ongrid, phi_ongrid, shift_ongrid)

    shift_ongrid = shift_ongrid .- mu

    # --------------- Causality --------------- #
    integrands = epsilon_helpers(electronic_spec, Z_ongrid, phi_ongrid, shift_ongrid, w_prime, idx_skip = [2])

    FLIP = ifelse.(integrands[1] .<= 0, 1, -1)

    # ------------- ε- & ω-integration ------------- #
    omega_A = integrands[3] .* FLIP .* tanhw     # ε + (χ(ω) - μ_F)  (shift_int branch)
    Ne_sc_A = dos_int + 2/π*trapz(w_prime, omega_A)
    root_eq = Ne_nsc - Ne_sc_A


    return root_eq
end


