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
function root_finding(fmu, outdir, fermi_level; shift = nothing, omega_shift = nothing, testMode::Bool = false)

    # Optional diagnostic dumped alongside muError when the update fails: the shift
    # channel χ(ω) over its frequency grid. N_e is obtained from an ω-integral of χ,
    # so if χ has not decayed to ~0 at the edge of its grid the integral is truncated
    # and the root find cannot converge (real axis: increase omega_c).
    function plotShift()
        testMode && return
        (isnothing(shift) || isnothing(omega_shift)) && return
        chi = real.(shift)
        plot(omega_shift, chi, label="χ(ω)", title="Shift channel used in the μ-update", xlabel="ω / meV", ylabel="χ / meV")
        savefig(outdir*"muError_shift.png")
        savePlotData(outdir*"muError_shift.dat", "#  ω / meV       χ(ω) / meV", omega_shift, chi)
    end

    ### starting values for mu
    mu_wndw = 50.0
    mu0 = fermi_level-mu_wndw
    mu1 = fermi_level+mu_wndw
    fmu0 = fmu(mu0)
    fmu1 = fmu(mu1)

    ### wrong slope 
    if fmu1 > fmu0
        # plot electron number vs mu
        mu_error = range(mu0 - mu_wndw, mu1 + mu_wndw, 50)
        Ne_error = zeros(size(mu_error))
        for k in eachindex(mu_error)
            Ne_error[k] = fmu(mu_error[k])
        end

        # diagnostics for the error below: not tied to flag_figure, but skipped in test
        # mode, where outdir is never created and savefig would mask the real message
        if ~testMode
            p = plot(mu_error, Ne_error, label="Ne_nsc - Ne_sc", title="Ne in normal state minus sc state")
            #vline(p, [mu0, mu1], label="mu")
            savefig(outdir*"muError.png")
            savePlotData(outdir*"muError.dat", "#  μ / meV       Ne_nsc - Ne_sc", mu_error, Ne_error)
        end
        plotShift()

        error("The number of electrons decreases with increasing mu! See muError.png/muError.dat, muError_shift.png (does χ(ω) decay to 0?) and the μ-update section of the Troubleshooting page.")
    end

    ### find minimum interval around ef in which a sign change occurs
    iter = 0
    mu_error = [mu0, mu1]
    fmu_error = [fmu0, fmu1]
    # the window is 2*mu_wndw wide, so it has to slide by its own width for the carried
    # over endpoint value to still belong to the μ it is assigned to: sliding by mu_wndw
    # would leave fmu1 holding f(mu0_old) while mu1 sits half a window further out
    mu_step = 2 * mu_wndw
    while fmu0 * fmu1 > 0
        if sign(fmu0) < 0
            mu0 -= mu_step
            mu1 -= mu_step
            fmu1 = fmu0
            fmu0 = fmu(mu0)
            pushfirst!(mu_error, mu0)
            pushfirst!(fmu_error, fmu0)
        else
            mu0 += mu_step
            mu1 += mu_step
            fmu0 = fmu1
            fmu1 = fmu(mu1)
            push!(mu_error, mu1)
            push!(fmu_error, fmu1)
        end
        iter += 1
        if iter > 50    # 50 * mu_step = 5 eV to either side
            if ~testMode
                plot(mu_error, fmu_error, label="Ne_nsc - Ne_sc", title="Ne in normal state minus sc state")
                savefig(outdir*"muError.png")
                savePlotData(outdir*"muError.dat", "#  μ / meV       Ne_nsc - Ne_sc", mu_error, fmu_error)
            end
            plotShift()

            mu0error = mu_error[1]
            mu1error = mu_error[end]
            error("Error in mu update - Couldn't find a root in the interval [$mu0error,$mu1error]. See muError.png/muError.dat, muError_shift.png (does χ(ω) decay to 0?) and the μ-update section of the Troubleshooting page.")
        end
    end

    ### calc new mu using the RegulaFalsi method
    mu = RegulaFalsi(fmu, fmu0, fmu1, mu0, mu1, 1e-3, 1e-6)     
    # Did some tests: ftol = 1e-4 would probably also be sufficient as mu only changes a few percent
    #                 This should not affect the Tc

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
function diff_Ne(mu, Ne_nsc, itemp, wsi, dos_en, dos, znormip, phiip, shiftip)
    # eq (9) & (11) overleaf
    diff = dos_en .- mu
    theta = (wsi' .* znormip') .^ 2 .+ (diff .+ shiftip') .^ 2 .+ phiip .^ 2
    summand_sc = (diff .+ shiftip') ./ theta .- diff ./ (wsi'.^ 2 .+ diff.^ 2)
    
    # matsubara sum
    summand_sc = dropdims((sum(summand_sc, dims=2)), dims=2)

    summand_sc = 2 * fermiFcn(dos_en, mu, itemp) - 4 * kb * itemp * summand_sc

    # diff between Ne in normal and sc state
    Ne_sc =  Ne_nsc - trapz(dos_en, dos.*summand_sc) 

    return Ne_sc

end


"""
    build_fmu_matsubara(itemp, wsi, dos_en, dos, znormip, phiphip, phicip, shiftip) -> fmu

Return the charge-neutrality residual `fmu(μ) = Nₑ_nsc - Nₑ_sc(μ)` used by the imaginary-axis
μ-update, *without* running the root find. This is the exact function whose root
[`update_mu_own`](@ref) solves; it is factored out so it can be scanned/plotted on its own for
debugging (see the μ-update testing notebook).
"""
function build_fmu_matsubara(itemp, wsi, dos_en, dos, znormip, phiphip, phicip, shiftip)

    # φ(ε, iωₙ) = φ_ph(iωₙ) + φ_c(ε) over the ε-window [-encut, encut]; φ_c is empty outside
    # vDOS+W, where φ carries no ε-dependence and a single row is enough. Bound to a new name
    # rather than reassigned: it is captured by `fmu` below, and assigning to a captured variable
    # boxes it, which would leave the closure - and every arithmetic on its result in the
    # root find - untyped. permutedims (not ') so that both branches give a Matrix{Float64}
    phi_cut = isempty(phicip) ? permutedims(phiphip) :
              permutedims(phiphip) .+ phicip

    ### Calculate N_e in the non-SC state
    Ne_nsc = trapz(dos_en, 2 .* fermiFcn(dos_en, 0.0, itemp) .* dos)

    # call calc_Ne_Sc with first argument unspecified
    fmu(x) = diff_Ne(x, Ne_nsc, itemp, wsi, dos_en, dos, znormip, phi_cut, shiftip)

    return fmu
end


"""
    update_mu_own(itemp, wsi, ef, dos_en, dos, znormip, phiphip, phicip, shiftip)

Routine to update chemical potential to fix the number of electrons
"""
function update_mu_own(itemp, wsi, dos_en, dos, znormip, phiphip, phicip, shiftip, fermi_level, outdir; testMode::Bool = false)

    fmu = build_fmu_matsubara(itemp, wsi, dos_en, dos, znormip, phiphip, phicip, shiftip)

    mu = root_finding(fmu, outdir, fermi_level; shift = shiftip, omega_shift = wsi, testMode = testMode)

    return mu

end



############################################################
# ---------------------- Real Axis ----------------------- #
############################################################
"""
    build_fmu_real_axis(itemp, w_static, w_prime, electronic_spec, dos_en, dos, znormip, phi_ph, phi_c, shiftip) -> fmu

Real-axis analogue of [`build_fmu_matsubara`](@ref): return the charge-neutrality residual
`fmu(μ) = Nₑ_nsc - Nₑ_sc(μ)` without running the root find, so it can be scanned/plotted for
debugging. The self-energy (`znormip`, `phi_ph`, `phi_c`, `shiftip`) is interpolated onto the
ω'-integration grid `w_prime`, so varying `w_prime` (or truncating `shiftip`/`w_static`) shows
directly how the μ-update integral reacts to the ω-grid.
"""
function build_fmu_real_axis(itemp, w_static, w_prime, electronic_spec, dos_en, dos, znormip, phi_ph, phi_c, shiftip)

    ### Calculate N_e in the non-SC state
    Ne_nsc = 2 .* trapz(dos_en, fermiFcn(dos_en, 0.0, itemp) .* dos)

    # interpolate Z, χ, φ_ph onto the ω'-grid (φ_c lives on the ε-grid, no ω-interp)
    Z_itp = linear_interpolation(w_static, znormip, extrapolation_bc=Flat())
    shift_itp = linear_interpolation(w_static, shiftip, extrapolation_bc=Flat())

    Z_ongrid = Z_itp.(w_prime)
    phi_ph_ongrid = linear_interpolation(w_static, phi_ph, extrapolation_bc=Flat()).(w_prime)
    shift_ongrid = shift_itp.(w_prime)

    dos_int = trapz(dos_en, dos)

    # temperature
    β = 1/(kb*itemp)
    tanhw = tanh.(β .* w_prime ./ 2)

    # call calc_Ne_Sc with first argument unspecified
    fmu(x)  = diff_Ne_realAxis(x, Ne_nsc, tanhw, w_prime, electronic_spec, dos_int, Z_ongrid, phi_ph_ongrid, phi_c, shift_ongrid)

    return fmu
end


"""
    mu_update_real_axis()

Update the chemical potential to conserve charge neutrality - real axis implementation.
φ is passed in separable form: φ_ph(ω) always, φ_c(ε) only in vDOS+W (empty in vDOS+μ).
"""
function mu_update_real_axis(itemp, fermi_level, w_static, w_prime, electronic_spec, dos_en, dos, znormip, phi_ph, phi_c, shiftip, outdir; testMode::Bool = false)

    fmu = build_fmu_real_axis(itemp, w_static, w_prime, electronic_spec, dos_en, dos, znormip, phi_ph, phi_c, shiftip)

    mu = root_finding(fmu, outdir, fermi_level; shift = shiftip, omega_shift = w_static, testMode = testMode)
    #@time mu2 = find_zero((fmu, dfmu), fermi_level, Roots.LithBoonkkampIJzerman(3, 1))

    return mu
end


"""
    diff_Ne_realAxis(mu, Ne_nsc, itemp, w_prime, dos_en, dos, znormip, deltaip, shiftip)

Calculate the number of electrons in the sc state for a given chemical
potential minus the number of electrons in the normal state
"""
function diff_Ne_realAxis(mu, Ne_nsc, tanhw, w_prime, electronic_spec, dos_int, Z_ongrid, phi_ph_ongrid, phi_c, shift_ongrid)

    shift_ongrid = shift_ongrid .- mu

    # Causality integrand. vDOS+μ (phi_c empty): φ(ω) ε-independent -> scalar epsilon_helpers.
    # vDOS+W: φ(ω,ε)=φ_ph(ω)+φ_c(ε) enters ε_p=√((ω'Z)²−φ²) per interval -> vDOSW moment loop
    # (Coulomb/W plays no role in charge conservation, so it is dropped here).
    if isempty(phi_c)
        integrands = epsilon_helpers(electronic_spec, Z_ongrid, phi_ph_ongrid, shift_ongrid, w_prime, idx_skip = [2])
        z_int, chi_int = integrands[1], integrands[3]
    else
        z_int, chi_int = epsilon_causality_vDOSW(electronic_spec, Z_ongrid, phi_ph_ongrid, phi_c, shift_ongrid, w_prime)
    end

    return Ne_root(Ne_nsc, dos_int, z_int, chi_int, tanhw, w_prime)
end

# Shared ε- & ω-integration: causality sign from the Z integrand, N_e from χ.
function Ne_root(Ne_nsc, dos_int, z_int, chi_int, tanhw, w_prime)
    FLIP = ifelse.(z_int .<= 0, 1, -1)
    omega_A = chi_int .* FLIP .* tanhw     # ε + (χ(ω) - μ_F)  (shift_int branch)
    Ne_sc_A = dos_int + 2/π*trapz(w_prime, omega_A)
    return Ne_nsc - Ne_sc_A
end


