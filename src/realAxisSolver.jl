"""
    File containing the real axis Eliashberg equations

"""


mutable struct RealAxisState
    Z::Vector{ComplexF64} 
    delta::Vector{ComplexF64}
    chi::Vector{ComplexF64}
    fermi_level::Float64
    phi_ph::Vector{ComplexF64}   # φ_ph(ω) on the w_static grid (both modes)
    phi_c::Vector{ComplexF64}    # φ_c(ε) on the dos_en grid (vDOS+W); empty in vDOS+μ
end

RealAxisState(Z, delta, chi, fermi_level) = RealAxisState(Z, delta, chi, fermi_level, ComplexF64[], ComplexF64[])


""" 

Real axis eliashberg solver.
"""
function RealAxisSolver(inp::arguments)
    # do not show plots
    default(show = false) 

    # placeholders: keep the error handler usable even if the run dies before the
    # log file exists
    log_file = IOBuffer()
    errorLogger = SimpleLogger(log_file, Logging.Error)
    Tc = [NaN, NaN]

    try
        dt = @elapsed begin

            strIsoME = printIsoME()

            ### Create directory
            inp, log_file, errorLogger = createDirectory(inp, strIsoME)

            ### Check input
            inp = stage("in input structure") do
                checkInput(inp, realSolver=true)
            end

            ### read inputs
            """
            ********** ToDo **********
            - adapt messages in printFlagsAsText
            - add names to start message
            """
            inp, console, matval, ML_Tc = stage("while reading the inputs") do
                InputParser(inp, log_file, mode=1)
            end

            ### Print to console ###
            printFlagsAsText(inp, log_file, mode="realFreq")

            ########### start loop over temperatures ##########
            Tc, temps, Znorm0, Delta0, Shift0 = stage("while solving the Eliashberg equations") do
                findTc_RealAxis(inp, console, matval, ML_Tc, log_file)
            end

            ### write Tc to console
            attempt(inp, log_file, "Error while printing the summary.") do
                printSummary(inp, Tc, log_file)
            end

            ### Outputs ###
            if ~inp.testMode # no output in test mode
                ### save inputs
                attempt(inp, log_file, "Error while creating the Info file.") do
                    createInfoFile(inp)
                end


                ### create Summary file
                attempt(inp, log_file, "Error while creating the Summary.dat file.") do
                    # sampled at the gap edge ω_g, not at ω = 0 - see `gap_edge_index`
                    header = "# T/K   Re{Δ(ω_g)}/meV   Im{Δ(ω_g)}/meV   Re{Z(ω_g)}/1   Im{Z(ω_g)}/1   "
                    out_vars = Array{Float64}(undef, length(Delta0), 5)
                    out_vars[:, 1] = temps
                    out_vars[:, 2] = real(Delta0)
                    out_vars[:, 3] = imag(Delta0)
                    out_vars[:, 4] = real(Znorm0)
                    out_vars[:, 5] = imag(Znorm0)
                    if inp.cDOS_flag == 0   # for later when vDOS is implemented
                        header = header * "Re{χ(ω_g)}/meV   Im{χ(ω_g)}/meV   ϵ_F-μ/meV   "
                        out_vars = hcat(out_vars, real(Shift0), imag(Shift0))
                    end
                    createSummaryFile(inp, Tc, out_vars, header)
                end


                """
                ---------------- ToDo -------------
                plot different figures, e.g. real and imag of Delta
                """
                ### figures
                if inp.flag_figure == 1
                    attempt(inp, log_file, "Error while plotting. Skipping plots.") do
                        createFigures(inp, matval, real(Delta0), temps, Tc, log_file)
                    end
                end
            end
        end


        # print time elapsed
        printTee(log_file, "\nTotal Runtime: " * string(round(dt, digits=2)) * " seconds\n")

    catch ex
        ### global error handler: the single place a fatal error is reported
        handleFatalError(ex, catch_backtrace(), inp, log_file, errorLogger)

        # the caller sees the original exception, not the IsoME wrapper
        throw(ex isa IsoMEError ? ex.cause : ex)
    finally
        # close & save - also on the way out of an error. sigint is disabled so
        # that an impatient second Ctrl+C can not interrupt the cleanup itself
        if ~inp.testMode
            Base.disable_sigint() do
                closeLog(log_file)
            end
        end
    end

    if inp.returnTc
        return Tc
    end

end


"""
    findTc(inp, console, matval, ML_Tc, log_file)

Solve the real axis eliashberg equations
"""
function findTc_RealAxis(inp, console, matval, ML_Tc, log_file)
    inp.temps = sort(inp.temps)
    nT = size(inp.temps, 1)
    Delta0 = Vector{ComplexF64}()
    Shift0 = Vector{ComplexF64}()
    Znorm0 = Vector{ComplexF64}()
    Tc = [NaN, NaN]

    realAxisState = nothing

    if inp.temps == [-1]    # Tc search mode
        # initial guess, Machine learning Tc
        itemp = maximum([1.0, round(ML_Tc)])

        # expansion of a + b*log(c-x) at x = 0, a=Delta(T2), b=1, c=Delta(T1)
        m(x, p) = p[1] + log(p[3]) .- p[2] * x ./ p[3] .- p[2] * x .^ 2 / (2 * p[3]^2) .- p[2] * x .^ 3 / (3 * p[3]^3) .- p[2] * x .^ 4 / (4 * p[3]^4) .- p[2] * x .^ 5 / (5 * p[3]^5)
        inp.temps = Vector{Float64}()
        fitFlag = true
        while true
            # save iterations
            inp.temps = push!(inp.temps, itemp)

            β = 1 / (kb * itemp)

            realAxisParameter = precompute(β, inp, matval, console, log_file)

            if inp.cDOS_flag == 0 && isnothing(realAxisState)
                # initial values vDOS
                realAxisState = initialize_real_axis_vDOS(itemp, inp, console.cDOS, matval, realAxisParameter, log_file)
            end

            # solve Eliashberg equations
            if inp.cDOS_flag == 0
                data, realAxisState = solve_realAxis_vDOS(itemp, inp, console.vDOS, matval, realAxisParameter, realAxisState, log_file)
            elseif inp.cDOS_flag == 1
                data, realAxisState = solve_realAxis_cDOS(itemp, inp, console.cDOS, matval, realAxisParameter, log_file)
            end
            if inp.cDOS_flag == 0
                Znorm0 = push!(Znorm0, data[1])
                Delta0 = push!(Delta0, data[2])
                Shift0 = push!(Shift0, data[3])
            elseif inp.cDOS_flag == 1
                Znorm0 = push!(Znorm0, data[1])
                Delta0 = push!(Delta0, data[2])
            end


            # Escape
            if itemp < 1 && isnan(Delta0[end])
                # log file
                text = "Lowest temperature of Tc search mode reached. To search at lower temperatures, set them manually.\n"
                print(log_file, text)
                print(text)

                Tc = [NaN, 0.5]
                break
            elseif length(inp.temps) > 500
                # log file
                print(log_file, "Couldn't find a Tc! \n")

                print("Couldn't find a Tc! \n")

                break
            end

            order = sortperm(inp.temps)
            if length(Delta0) >= 2 && any(diff(inp.temps[order]) .<= 1 .& (.~isnan.(Delta0[order][1:end-1]) .& isnan.(Delta0[order][2:end])))
                # converged and sort
                inp.temps = inp.temps[order]
                Delta0 = Delta0[order]
                Tc = [maximum(inp.temps[.~isnan.(Delta0)]), minimum(inp.temps[isnan.(Delta0)])]
                break

            elseif all(isnan.(Delta0))
                # temperature too high
                if itemp > 1
                    itemp = ceil(itemp / 2)
                else
                    # lowest T
                    itemp = 1 / 2
                end

            elseif sum(.~isnan.(Delta0)) == 1
                # get a second gap value
                if length(Delta0) == 1
                    itemp = itemp + round(maximum([itemp / 2, 2]))
                else
                    itemp = ceil((maximum(inp.temps[.~isnan.(Delta0)]) + minimum(inp.temps[isnan.(Delta0)])) / 2)
                end

            elseif sum(.~isnan.(Delta0)) == 2 && fitFlag
                # fit only once
                fitFlag = false

                nnanDelta = .~isnan.(Delta0)
                realDelta0 = real(Delta0[nnanDelta])
                # fit gap values
                p0 = convert(Vector{Float64}, [maximum(realDelta0), 1, minimum(realDelta0)])
                try
                    fit = curve_fit(m, inp.temps[nnanDelta], realDelta0, p0)
                    par = fit.param

                    # find root
                    m2(x) = m(x, par)
                    itemp = floor(find_zero(m2, par[2]))

                catch ex
                    ex isa InterruptException && rethrow(ex)
                    # expansion to third order --> analytical formula for root (only one real root)
                    a = p0[1]
                    c = p0[3]
                    itemp = -c / 2
                    itemp += -(3 * c^2) / (2 * (5 * c^3 + 12 * a * c^3 +
                                                2 * sqrt(13 * c^6 + 30 * a * c^6 + 36 * a^2 * c^6))^(1 / 3))
                    itemp += 0.5 * (5 * c^3 + 12 * a * c^3 +
                                    2 * sqrt(13 * c^6 + 30 * a * c^6 + 36 * a^2 * c^6))^(1 / 3)
                    itemp = round(itemp)
                    # formula works only if a,c are far away from the Tc
                    # if a,c are close to Tc the estimated T will be too small but this case is caputred by the sanity check
                end

                # sanity check
                if length(Delta0) == 2  # no nans
                    if itemp < maximum(inp.temps)       # error in fit
                        itemp = ceil(maximum(inp.temps) * 3 / 2)
                    elseif itemp == maximum(inp.temps)  # fit equals highest converged value
                        itemp += maximum([2, round(itemp / 10)])
                    end
                elseif length(Delta0) >= 2
                    # prevent search below converged T, above not converged T
                    if itemp <= maximum(inp.temps[nnanDelta]) || itemp >= minimum(inp.temps[isnan.(Delta0)])
                        itemp = ceil((maximum(inp.temps[nnanDelta]) + minimum(inp.temps[isnan.(Delta0)])) / 2)
                    end

                end
            else
                # search around fit value
                if any(isnan.(Delta0))
                    itemp = ceil((maximum(inp.temps[.~isnan.(Delta0)]) + minimum(inp.temps[isnan.(Delta0)])) / 2)
                else
                    itemp += maximum([2, round(itemp / 5)])
                end
            end

        end

    else    # given temperatures
        for iT in 1:nT
            itemp = inp.temps[iT]

            β = 1 / (kb * itemp)

            realAxisParameter = precompute(β, inp, matval, console, log_file)


            # initial values vDOS
            if inp.cDOS_flag == 0 && iT == 1
                realAxisState = initialize_real_axis_vDOS(itemp, inp, console.cDOS, matval, realAxisParameter, log_file)
            end

            # solve Eliashberg equations
            if inp.cDOS_flag == 0
                data, realAxisState = solve_realAxis_vDOS(itemp, inp, console.vDOS, matval, realAxisParameter, realAxisState, log_file)
            elseif inp.cDOS_flag == 1
                data, realAxisState = solve_realAxis_cDOS(itemp, inp, console.cDOS, matval, realAxisParameter, log_file)
            end
            if inp.cDOS_flag == 0
                Znorm0 = push!(Znorm0, data[1])
                Delta0 = push!(Delta0, data[2])
                Shift0 = push!(Shift0, data[3])
            elseif inp.cDOS_flag == 1
                Znorm0 = push!(Znorm0, data[1])
                Delta0 = push!(Delta0, data[2])
            end

            if isnan(Delta0[end])
                # escape
                inp.temps = inp.temps[1:length(Delta0)]
                if length(inp.temps) > 1
                    Tc[1] = inp.temps[end-1]
                end
                # upper bound Tc
                Tc[2] = inp.temps[end]
                break
            end
        end

        if all(.~isnan.(Delta0))
            # lower bound Tc
            Tc[1] = inp.temps[end]
        end

    end

    printTextCentered("Stopping now!", console.partingLine, file=log_file, bold=true)

    return Tc, inp.temps, Znorm0, Delta0, Shift0

end



"""
    solve_realAxis_cDOS(itemp, inp, console, realAxisParameter, log_file; vDOS_initial_guess=false)

Solve the real-axis Eliashberg equations in the cDOS+μ approximation.
When `vDOS_initial_guess=true`, the cDOS solution is returned in a `RealAxisState`
with zero χ so it can seed the vDOS solver.
"""
function solve_realAxis_cDOS(itemp, inp, console, matval, realAxisParameter, log_file; vDOS_initial_guess::Bool=false)
    (; reOmega_c, N_it, conv_thr, minGap, nItFullCoul, min_it) = inp
    (w_static, _) = realAxisParameter
    (_, _, _, _, Weep, dosef, idx_ef, _, BCS_gap, _) = matval
    # When used as vDOS+W initializer, muc_ME is unset (nothing). Derive it from W(εF,εF)·N(εF)
    # using the Morel-Anderson formula, consistent with calcMucs() in ReadIn.jl.
    muc_ME = if vDOS_initial_guess && isnan(inp.muc_ME) && inp.include_Weep == 1
        typEl = !isnan(inp.typEl) ? inp.typEl : inp.efW
        mu = Weep[idx_ef, idx_ef] * dosef
        mu / (1 + mu * log(typEl / reOmega_c))
    else
        Float64(inp.muc_ME)
    end

    # smaller threshold for initial guess calculation
    conv_thr = vDOS_initial_guess ? max(1e-3, conv_thr) : conv_thr

    title = vDOS_initial_guess ? "Initial guess: cDOS at T = " * string(itemp) * " K" : "T = " * string(itemp) * " K "
    printTextCentered(title, console.partingLine, file=log_file, bold=true)
    printTee(log_file, "\n")

    state = initial_real_axis_state(inp, realAxisParameter, BCS_gap)
    ig0 = gap_edge_index(w_static, BCS_gap)   # same sampling point as the iteration rows
    initValues = [0 real(state.Z[ig0]) imag(state.Z[ig0]) real(state.delta[ig0]) imag(state.delta[ig0]) nothing]
    printTableHeader(console, initValues, log_file)

    β = 1 / (kb * itemp)
    Z_new = state.Z
    delta_new = state.delta
    data = [Z_new[ig0], delta_new[ig0]]

    # --- head/tail (vDOS-style) comparison solver ---
    # pole-anchored Chebyshev head + static linear tail; the tail kernels are precomputed
    # once and the head is rebuilt only when the pole of Θ(ω') moves (maybe_refresh_head!).
    ws_ht = build_cDOS_wprime_workspace(inp, realAxisParameter, BCS_gap)

    # seed for the gap estimate; from here on it is the root of Θ(ω) = ω − ReΔ(ω),
    # never Δ(w_static[1]) - see spectral_gap in wprimeGrid.jl
    gap0 = BCS_gap

    # Iterate
    for i_it in 1:N_it
        delta_prev = copy(delta_new)
        Z_prev = copy(Z_new)
        broyden_beta = mixing_parameter(inp, i_it)

        gap0 = spectral_gap(w_static, delta_prev, gap0)
        pole = gap0
        if pole >= ws_ht.wp_max
            ws_ht = build_cDOS_wprime_workspace(inp, realAxisParameter, pole)
        end
        maybe_refresh_head!(ws_ht, [pole], w_static, inp.n_cheb)

        # weight coulomb interaction
        wgCoulomb = minimum([1, i_it / nItFullCoul])

        # Update Eliashberg
        Z_new, delta_new = realEliashbergEq(muc_ME, β, delta_prev, ws_ht, w_static, wgCoulomb)

        Z_new = (1.0 - abs(broyden_beta)) .* Z_prev .+ abs(broyden_beta) .* Z_new
        delta_new = (1.0 - abs(broyden_beta)) .* delta_prev .+ abs(broyden_beta) .* delta_new


        rel_delta = sum(abs.(delta_new - delta_prev))
        abs_delta = sum(abs.(delta_new))
        convergence = rel_delta / abs_delta
        #convergence = abs(sqrt(sum(abs2.(delta_new .- delta_prev))/length(delta_new))/gap0)  # Alejandros criteria
        # Z and Δ are taken at the gap edge, not at w_static[1] - see `gap_edge_index`
        idx_gapEdge = gap_edge_index(w_static, gap0)
        data = [Z_new[idx_gapEdge], delta_new[idx_gapEdge]]
        outputVec = real_axis_cdOS_output(i_it, Z_new, delta_new, convergence, idx_gapEdge)
        nan_state = any(isnan, outputVec)
        print_real_axis_iteration(outputVec, console, log_file)

        if abs(convergence) < conv_thr && i_it > maximum([min_it, nItFullCoul + 1])
            print_real_axis_converged(itemp, console, log_file)

            state.Z = Z_new
            state.delta = delta_new
            state.phi_ph = delta_new .* Z_new
            if !vDOS_initial_guess
                save_real_axis_cDOS_outputs(itemp, inp, state, w_static, log_file)
            end
            return data, state
        end

        if real(data[2]) < minGap && i_it > maximum([min_it, nItFullCoul + 1])
            print_real_axis_gap_too_small(itemp, minGap, console, log_file)
            data[2] = NaN

            # initial guess vDOS
            if  vDOS_initial_guess
                state.Z = Z_new
                state.delta = delta_new
                state.phi_ph = delta_new .* Z_new
            end

            return data, state
        end

        # A NaN anywhere in the console row abandons this temperature. Every self-energy
        # component is built from sums over the whole ω'-grid, so a single NaN entry
        # contaminates all of them within one iteration and iterating on is pointless.
        # Handled like the N_it case (Δ = NaN), so the Tc search continues at lower T.
        if i_it == N_it || nan_state
            nan_state ? print_real_axis_nan(itemp, console, log_file) :
                        print_real_axis_not_converged(inp, console, log_file)
            data[2] = NaN

            # initial guess vDOS; a NaN state must not be handed on - the vDOS solver
            # would never recover from it, so leave the last finite guess in place.
            if vDOS_initial_guess && !nan_state
                state.Z = Z_new
                state.delta = delta_new
                state.phi_ph = delta_new .* Z_new
            end

            return data, state
        end
    end
end

"""
    solve_realAxis_vDOS(itemp, inp, console, matval, realAxisParameter, state, log_file)

Solve the real-axis Eliashberg equations in the vDOS+μ approximation.
"""
function solve_realAxis_vDOS(itemp, inp, console, matval, realAxisParameter, state::RealAxisState, log_file)
    # destruct inputs
    (_, _, dos_en, dos, Weep, dosef, idx_ef, _, BCS_gap, _) = matval
    (; muc_ME, mu_flag, N_it, conv_thr, minGap, nItFullCoul, min_it) = inp
    (w_static, _) = realAxisParameter


    printTextCentered("T = " * string(itemp) * " K ", console.partingLine, file=log_file, bold=true)
    printTee(log_file, "\n")

    # ----- electronic spectral information ----- #
    # INFO: Could be calculated in read in once instead of here
    eps_1, eps_2 = transpose(dos_en[1:end-1]), transpose(dos_en[2:end])
    dos_1, dos_2 = transpose(dos[1:end-1]), transpose(dos[2:end])
    dε = eps_2 .- eps_1
    ddos = dos_2 .- dos_1

    if inp.include_Weep == 1
        # Coulomb spectral integral: piecewise-linear moments of N(ε')W(ε,ε')
        NW0 = Weep[:, 1:end-1] .* dos_1
        NW1 = Weep[:, 2:end] .* dos_2
        M1W = (NW1 .- NW0) ./ dε
        M0W = NW0 .- eps_1 .* M1W
    else
        M1W = nothing
        M0W = nothing
    end

    electronic_spec = (eps_1, eps_2, dos_1, dos_2, dε, ddos, M0W, M1W)

    # ----- initial guesses -----#
    Z_new = state.Z
    delta_new = state.delta
    chi_new = state.chi
    # φ is the primary sc variable handed to realEliashbergEq; Δ is derived as φ/Z after mixing.
    # vDOS+W: φ(ω,ε) is an nω×ndos matrix; vDOS+μ: φ(ω) = Δ(ω)·Z(ω) is a vector.
    # φ_ph(ω): previous state, or the cDOS initial guess φ = Δ·Z
    phi_ph_new = isempty(state.phi_ph) ? (delta_new .* Z_new) : state.phi_ph
    # φ_c(ε): ε-resolved Coulomb, vDOS+W only (empty in vDOS+μ); zero on the first guess
    phi_c_new = if inp.include_Weep == 1
        isempty(state.phi_c) ? zeros(ComplexF64, length(dos_en)) : state.phi_c
    else
        ComplexF64[]
    end
    fermi_level = state.fermi_level

    # print console table; sampled at the gap edge, like the iteration rows
    ig0 = gap_edge_index(w_static, spectral_gap(w_static, state.delta, BCS_gap))
    initValues = [0 real(state.Z[ig0]) imag(state.Z[ig0]) real(state.chi[ig0]) imag(state.chi[ig0]) state.fermi_level real(state.delta[ig0]) imag(state.delta[ig0]) nothing]
    printTableHeader(console, initValues, log_file)

    β = 1 / (kb * itemp)
    data = [state.Z[ig0], state.delta[ig0], state.chi[ig0]]

    # ------ integration grid and kernels ----- #
    # head/tail split (vDOS): tail kernels are computed once per temperature,
    # the pole-anchored head grid is rebuilt only when the poles move.
    # wp_max = 2·Δ(0) of the starting state (cDOS solution / previous temperature).
    gridws = build_wprime_workspace(inp, realAxisParameter, spectral_gap(w_static, state.delta, BCS_gap))

    # seed for the gap estimate; from here on it is the root of Θ(ω) = ω − ReΔ(ω),
    # never Δ(w_static[1]) - see spectral_gap in wprimeGrid.jl
    gap0 = spectral_gap(w_static, state.delta, BCS_gap)

    # optional Broyden mixer (broyden_flag == 1); default is linear mixing
    broyden = inp.broyden_flag == 1 ? BroydenMixer(inp.broyden_mem) : nothing

    for i_it in 1:N_it
        delta_prev = copy(delta_new)
        Z_prev = copy(Z_new)
        chi_prev = copy(chi_new)
        phi_ph_prev = copy(phi_ph_new)
        phi_c_prev = copy(phi_c_new)
        broyden_beta = mixing_parameter(inp, i_it)
        gap0 = spectral_gap(w_static, delta_prev, gap0)

        # locate the ω'-integrand poles (vDOS+W: from the modified S/P quantities)
        if inp.include_Weep == 1
            poles = find_integrand_poles_vDOSW(gridws.wp_full, w_static, Z_prev, phi_ph_prev, phi_c_prev, chi_prev, fermi_level, dos_en, idx_ef, gap0)
        else
            poles = find_integrand_poles(gridws.wp_full, w_static, Z_prev, delta_prev, chi_prev, fermi_level, gap0)
        end
        if maximum(poles) >= gridws.wp_max
            # a pole moved past the fixed tail: rebuild with a larger head region
            gridws = build_wprime_workspace(inp, realAxisParameter, maximum(poles))
        end
        maybe_refresh_head!(gridws, poles, w_static, inp.n_cheb)
        w_prime = gridws.wp_full

        # weight coulomb interaction
        wgCoulomb = minimum([1, i_it / nItFullCoul])


        # Add mu thermalization: start mu-update only after a few iterations
        if mu_flag == 1 # && i_it > maximum([min_it, nItFullCoul + 1]) -1
            fermi_level = mu_update_real_axis(itemp, fermi_level, w_static, w_prime, electronic_spec, dos_en, dos, Z_prev, phi_ph_prev, phi_c_prev, chi_prev, inp.outdir)
        end

        if inp.include_Weep == 1
            Z_new, chi_new, phi_ph_new, phi_c_new = realEliashbergEq(β, Z_prev, phi_ph_prev, phi_c_prev, chi_prev, gridws, w_static, dosef, fermi_level, wgCoulomb, electronic_spec)
        else
            Z_new, chi_new, phi_ph_new = realEliashbergEq(muc_ME, β, Z_prev, phi_ph_prev, chi_prev, gridws, w_static, dosef, wgCoulomb, fermi_level, electronic_spec)
            phi_c_new = ComplexF64[]
        end

        # mixing on the primary variables (Z, χ, φ_ph, φ_c); φ_c is empty in vDOS+μ (no-op)
        if inp.broyden_flag == 1
            Z_new, chi_new, phi_ph_new, phi_c_new = broyden_mix!(broyden, abs(broyden_beta), Z_prev, chi_prev, phi_ph_prev, phi_c_prev, Z_new, chi_new, phi_ph_new, phi_c_new)
        else
            β_mix = abs(broyden_beta)
            chi_new    = (1.0 - β_mix) .* chi_prev    .+ β_mix .* chi_new
            Z_new      = (1.0 - β_mix) .* Z_prev      .+ β_mix .* Z_new
            phi_ph_new = (1.0 - β_mix) .* phi_ph_prev .+ β_mix .* phi_ph_new
            phi_c_new  = (1.0 - β_mix) .* phi_c_prev  .+ β_mix .* phi_c_new
        end

        # Δ(ω) = φ(ω, ε_F)/Z(ω); φ(ω, ε_F) = φ_ph(ω) + φ_c(ε_F) (vDOS+W) or φ_ph(ω) (vDOS+μ)
        phi_ef = inp.include_Weep == 1 ? phi_ph_new .+ phi_c_new[idx_ef] : phi_ph_new
        delta_new = phi_ef ./ Z_new

        convergence = sqrt(sum(abs2.(delta_new .- delta_prev))/length(delta_new))
        # Z, Δ and χ are taken at the gap edge, not at w_static[1] - see `gap_edge_index`
        idx_gapEdge = gap_edge_index(w_static, gap0)
        data = [Z_new[idx_gapEdge], delta_new[idx_gapEdge], chi_new[idx_gapEdge]]
        outputVec = real_axis_vDOS_output(i_it, Z_new, delta_new, chi_new, fermi_level, convergence, idx_gapEdge)
        nan_state = any(isnan, outputVec)
        print_real_axis_iteration(outputVec, console, log_file)

        if abs(convergence) < conv_thr && i_it > maximum([min_it, nItFullCoul + 1])
            print_real_axis_converged(itemp, console, log_file)

            if inp.flag_writeSelfEnergy == 1
                try
                    saveSelfEnergy_realAxis(itemp, inp, collect(w_static), delta_new, Z_new,
                                            chi=chi_new,
                                            phi_ph=(inp.include_Weep == 1 ? phi_ph_new : nothing),
                                            phi_c=(inp.include_Weep == 1 ? phi_c_new : nothing),
                                            epsilon=(inp.include_Weep == 1 ? dos_en : nothing))
                catch ex
                    writeToCrashFile(inp)
                    printWarning("Error while saving self energy components.", log_file, ex=ex)
                end
            end

            # next guess
            # state.Z = Z_new
            # state.chi = chi_new
            # state.delta = delta_new
            # state.phi_ph = phi_ph_new
            # state.phi_c = phi_c_new
            # state.fermi_level = fermi_level
            return data, state
        end

        if real(data[2]) < minGap && i_it > maximum([min_it, nItFullCoul + 1])
            print_real_axis_gap_too_small(itemp, minGap, console, log_file)
            data[2] = NaN
            return data, state
        end

        # A NaN anywhere in the console row abandons this temperature - see the cDOS loop.
        # `state` is deliberately left untouched, so the next temperature restarts from the
        # last converged self-energy rather than from the NaN one.
        if i_it == N_it || nan_state
            nan_state ? print_real_axis_nan(itemp, console, log_file) :
                        print_real_axis_not_converged(inp, console, log_file)
            data[2] = NaN
            return data, state
        end
    end
end


##############################################################
# -------------------- Helper functions -------------------- #
##############################################################
function initial_real_axis_state(inp, realAxisParameter, BCS_gap)
    (w_static, _) = realAxisParameter
    nw = length(w_static)

    return RealAxisState(
        ones(ComplexF64, nw),                        # Z
        ones(ComplexF64, nw) .* (BCS_gap + im * 1e-4),   # Delta
        -zeros(ComplexF64, length(w_static)),               # Chi
        0.0,                                                # fermi-level
        ComplexF64[],                                       # φ_ph(ω)
        ComplexF64[],                                       # φ_c(ε), vDOS+W only
    )
end

"""
    mixing_parameter(inp, i_it)

linear mixing factor of eliashberg solutions
"""
function mixing_parameter(inp, i_it)
    if isnan(inp.mixing_beta)
        return maximum([0.5, 1.0 - 0.05 * (i_it - 1)])
    end

    return inp.mixing_beta
end


function save_real_axis_cDOS_outputs(itemp, inp, state, w_static, log_file)
    if inp.flag_writeSelfEnergy == 1
        try
            saveSelfEnergy_realAxis(itemp, inp, collect(w_static), state.delta, state.Z)
        catch ex
            writeToCrashFile(inp)
            printWarning("Error while saving self energy components.", log_file, ex=ex)
        end
    end
end


"""
    initialize_real_axis_vDOS(itemp, inp, console, realAxisParameter, log_file)

Perform a cDOS calculation as initial guess for vDOS.
"""
function initialize_real_axis_vDOS(itemp, inp, console, matval, realAxisParameter, log_file)
    _, state = solve_realAxis_cDOS(itemp, inp, console, matval, realAxisParameter, log_file; vDOS_initial_guess=true)
    return state
end




#################################################################
# ----------------------- print Helpers ----------------------- #
#################################################################
"""
    print_real_axis_iteration(outputVec, console, log_file)

print state of current iteration
"""
function print_real_axis_iteration(outputVec, console, log_file)
    outputVec, strConsole, format = formatTableRow(outputVec, console.width, console.precision)

    printTableRow(stdout, outputVec, strConsole, format)
    printTableRow(log_file, outputVec, strConsole, format)
end


"""
    print_real_axis_converged(itemp, console, log_file)

Print temperature converged
"""
function print_real_axis_converged(itemp, console, log_file)
    println(replace(console.Hline, "." => " "))
    printstyled("\nConvergence achieved for T = " * string(itemp) * " K\n"; bold=false)

    println(log_file, replace(console.Hline, "." => " "))
    printstyled(log_file, "\nConvergence achieved for T = " * string(itemp) * " K\n"; bold=false)
end


"""
    print_real_axis_gap_too_small(itemp, minGap, console, log_file)

gap at temperature too small
"""
function print_real_axis_gap_too_small(itemp, minGap, console, log_file)
    println(replace(console.Hline, "." => " "))
    printstyled("\nTemperature (T = " * string(itemp) * " K) too high, gap value already smaller than " * string(round(minGap, digits=2)) * " meV!\n\n"; bold=false)

    println(log_file, replace(console.Hline, "." => " "))
    printstyled(log_file, "\nTemperature (T = " * string(itemp) * " K) too high, gap value already smaller than " * string(round(minGap, digits=2)) * " meV!\n\n"; bold=false)
end


"""
    print_real_axis_nan(itemp, console, log_file)

self energy turned NaN, temperature abandoned
"""
function print_real_axis_nan(itemp, console, log_file)
    println(replace(console.Hline, "." => " "))
    printstyled("\nSelf energy became NaN at T = " * string(itemp) * " K!\n\n"; bold=true)

    println(log_file, replace(console.Hline, "." => " "))
    printstyled(log_file, "\nSelf energy became NaN at T = " * string(itemp) * " K!\n\n"; bold=true)
end


"""
     print_real_axis_not_converged(inp, console, log_file)

print max iterations exceeded
"""
function print_real_axis_not_converged(inp, console, log_file)
    println(replace(console.Hline, "." => " "))
    printstyled("\nConvergence not achieved within " * string(inp.N_it) * " iterations\n"; bold=true)
    println("\n")

    println(log_file, replace(console.Hline, "." => " "))
    printstyled(log_file, "\nConvergence not achieved within " * string(inp.N_it) * " iterations\n"; bold=true)
    println(log_file, "\n")
end

"""
    real_axis_cdOS_output(i_it, Z_new, delta_new, convergence, ig)

cDOS console output. Z and Δ are reported at the gap-edge index `ig` (see `gap_edge_index`),
the same point the solver uses for its own gap check, and the error column is the quantity
that is tested against `conv_thr`.
"""
function real_axis_cdOS_output(i_it, Z_new, delta_new, convergence, ig)
    return [i_it, real(Z_new[ig]), imag(Z_new[ig]), real(delta_new[ig]), imag(delta_new[ig]), abs(convergence)]
end


"""
    real_axis_vDOS_output(i_it, Z_new, delta_new, chi_new, fermi_level, convergence, ig)

vDOS console output. Z, χ and Δ are reported at the gap-edge index `ig` (see `gap_edge_index`),
the same point the solver uses for its own gap check, and the error column is the quantity
that is tested against `conv_thr`.
"""
function real_axis_vDOS_output(i_it, Z_new, delta_new, chi_new, fermi_level, convergence, ig)
    return [i_it, real(Z_new[ig]), imag(Z_new[ig]), real(chi_new[ig]), imag(chi_new[ig]), fermi_level,
            real(delta_new[ig]), imag(delta_new[ig]), abs(convergence)]
end


"""
    gap_edge_index(w_static, gap0)

Index of the ω-grid point closest to the gap edge `gap0`, i.e. where the console row is
sampled. Falls back to the first point if `gap0` is not usable.

The row must not be read off at `w_static[1]`: at finite temperature Z(ω) = 1 − Iz(ω)/ω
diverges as 1/ω (the thermal quasiparticle scattering rate is nonzero), so Z is huge and
Δ = φ/Z is crushed at the first grid point - for H3S at 180 K the table showed
Im Z = 37.8 and Δ = 0.30 meV against a gap of ~66 meV. Both are artefacts of ω₁ = 0.1 meV
being an arbitrary grid choice, and both vanish at low temperature, which is why this only
ever looked wrong for high-Tc materials.
"""
function gap_edge_index(w_static, gap0::Float64)
    isfinite(gap0) || return firstindex(w_static)
    idx = firstindex(w_static)
    best = Inf
    @inbounds for i in eachindex(w_static)
        d = abs(w_static[i] - gap0)
        if d < best
            best = d
            idx = i
        end
    end
    return idx
end





"""

Extract numbers from string
"""
function numbersFromString(str::String)
    nums = []
    current_number = ""

    for char in str
        if isdigit(char) || char == '.'  # Check if the character is a digit or decimal point
            current_number *= char      # Build the number string
        elseif current_number != ""     # If we encounter a non-digit and have a number built
            push!(nums, parse(Float64, current_number))  # Convert and store the number
            current_number = ""  # Reset for the next number
        end
    end

    # If a number is left at the end of the string
    if current_number != ""
        push!(nums, parse(Float64, current_number))
    end

    return Float64.(nums)
end
