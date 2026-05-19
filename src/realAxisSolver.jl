"""
    File containing the real axis Eliashberg equations

"""


mutable struct RealAxisState
    Z::Vector{ComplexF64} 
    delta::Vector{ComplexF64}
    chi::Vector{ComplexF64}
    fermi_level::Float64
    phi::Union{Nothing, Matrix{ComplexF64}}
end

RealAxisState(Z, delta, chi, fermi_level) = RealAxisState(Z, delta, chi, fermi_level, nothing)


""" 

Real axis eliashberg solver.
"""
function RealAxisSolver(inp::arguments)

    dt = @elapsed begin

        strIsoME = printIsoME()

        ### Create directory
        inp, log_file, errorLogger = createDirectory(inp, strIsoME)

        ### Check input
        try
            inp = checkInput(inp, realSolver=true)

        catch ex
            # crash file
            writeToCrashFile(inp)

            # console / log file
            printError("in input structure. Stopping now!", ex, log_file, errorLogger)

            rethrow(ex)
        end

        ### read inputs
        matval = ()
        ML_Tc = NaN
        console = Dict()
        a2F_itp = nothing
        try
            """ 
            ********** ToDo **********
            - adapt messages in printFlagsAsText
            - add names to start message
            """
            inp, console, matval, ML_Tc, a2F_itp = InputParser(inp, log_file, mode=1)
        catch ex
            # crash file
            writeToCrashFile(inp)

            # console / log file
            printError("while reading the inputs. Stopping now!", ex, log_file, errorLogger)

            rethrow(ex)
        end

        ### Print to console ###
        printFlagsAsText(inp, log_file, mode="realFreq")

        ########### start loop over temperatures ##########
        Tc = [NaN, NaN]
        temps = Vector{Float64}()
        Delta0 = Vector{Float64}()
        Shift0 = Vector{Float64}()
        Znorm0 = Vector{Float64}()

        try
            Tc, temps, Znorm0, Delta0, Shift0 = findTc_RealAxis(inp, console, matval, ML_Tc, a2F_itp, log_file)
        catch ex
            # crash file
            writeToCrashFile(inp)

            # console / log file
            printError("while solving the Eliashberg equations. Stopping now!", ex, log_file, errorLogger)
        end

        ### write Tc to console
        try
            printSummary(inp, Tc, log_file)
        catch ex
            # crash file
            writeToCrashFile(inp)

            # console / log file
            printWarning("Error while printing the summary.", log_file, ex=ex)
        end

        ### Outputs ###
        if ~inp.testMode # no output in test mode
            ### save inputs
            try
                createInfoFile(inp)
            catch ex
                # crash file
                writeToCrashFile(inp)

                # console / log file
                printWarning("Error while creating the Info file.", log_file, ex=ex)
            end


            ### create Summary file
            header = "# T/K   Re{Δ(0)}/meV   Im{Δ(0)}/meV   Re{Z(0)}/1   Im{Z(0)}/1   "
            out_vars = Array{Float64}(undef, length(Delta0), 5)
            out_vars[:, 1] = temps
            out_vars[:, 2] = real(Delta0)
            out_vars[:, 3] = imag(Delta0)
            out_vars[:, 4] = real(Znorm0)
            out_vars[:, 5] = imag(Znorm0)
            if inp.cDOS_flag == 0   # for later when vDOS is implemented
                header = header * "Re{χ(0)}/meV   Im{χ(0)}/meV   ϵ_F-μ/meV   "
                out_vars = hcat(out_vars, real(Shift0), imag(Shift0))
            end
            try
                createSummaryFile(inp, Tc, out_vars, header)
            catch ex
                # crash file
                writeToCrashFile(inp)

                # console / log file
                printWarning("Error while creating the Summary.dat file.", log_file, ex=ex)
            end


            """
            ---------------- ToDo -------------
            plot different figures, e.g. real and imag of Delta 
            """
            ### figures
            if inp.flag_figure == 1
                try
                    createFigures(inp, matval, real(Delta0), temps, Tc, log_file)
                catch ex
                    # crash file
                    writeToCrashFile(inp)

                    # console / log file
                    printWarning("Error while plotting. Skipping plots.", log_file, ex=ex)
                end
            end
        end
    end


    # print time elapsed
    print("\nTotal Runtime: ", round(dt, digits=2), " seconds\n")
    # log file
    print(log_file, "\nTotal Runtime: ", round(dt, digits=2), " seconds\n")

    # close & save
    if ~inp.testMode
        close(log_file)
    end

    if inp.returnTc
        return Tc
    end

end


"""
    findTc(inp, console, matval, ML_Tc, log_file)

Solve the real axis eliashberg equations
"""
function findTc_RealAxis(inp, console, matval, ML_Tc, a2F_itp, log_file)
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

            realAxisParameter = precompute(β, inp, matval, a2F_itp, console, log_file)

            if inp.cDOS_flag == 0 && isnothing(realAxisState)
                # initial values vDOS
                realAxisState = initialize_real_axis_vDOS(itemp, inp, console["cDOS"], realAxisParameter, log_file)
            end

            # solve Eliashberg equations
            if inp.cDOS_flag == 0
                data, realAxisState = solve_realAxis_vDOS(itemp, inp, console["vDOS"], matval, realAxisParameter, realAxisState, log_file)
            elseif inp.cDOS_flag == 1
                data, realAxisState = solve_realAxis_cDOS(itemp, inp, console["cDOS"], realAxisParameter, log_file)
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

                catch
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

            realAxisParameter = precompute(β, inp, matval, a2F_itp, console, log_file)


            # initial values vDOS
            if inp.cDOS_flag == 0 && iT == 1
                realAxisState = initialize_real_axis_vDOS(itemp, inp, console["cDOS"], realAxisParameter, log_file)
            end

            # solve Eliashberg equations
            if inp.cDOS_flag == 0
                data, realAxisState = solve_realAxis_vDOS(itemp, inp, console["vDOS"], matval, realAxisParameter, realAxisState, log_file)
            elseif inp.cDOS_flag == 1
                data, realAxisState = solve_realAxis_cDOS(itemp, inp, console["cDOS"], realAxisParameter, log_file)
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

    printTextCentered("Stopping now!", console["partingLine"], file=log_file, bold=true)

    return Tc, inp.temps, Znorm0, Delta0, Shift0

end



"""
    solve_realAxis_cDOS(itemp, inp, console, realAxisParameter, log_file; vDOS_initial_guess=false)

Solve the real-axis Eliashberg equations in the cDOS+μ approximation.
When `vDOS_initial_guess=true`, the cDOS solution is returned in a `RealAxisState`
with zero χ so it can seed the vDOS solver.
"""
function solve_realAxis_cDOS(itemp, inp, console, realAxisParameter, log_file; vDOS_initial_guess::Bool=false)
    (; reOmega_c, muc_ME, N_it, conv_thr, minGap, nItFullCoul, min_it) = inp
    (Kp_func, Km_func, w_static, w_dynam, _) = realAxisParameter

    # smaller threshold for initial guess calculation
    conv_thr = vDOS_initial_guess ? max(1e-3, conv_thr) : conv_thr

    title = vDOS_initial_guess ? "Initial guess: cDOS at T = " * string(itemp) * " K" : "T = " * string(itemp) * " K "
    printTextCentered(title, console["partingLine"], file=log_file, bold=true)
    printTee(log_file, "\n")

    state = initial_real_axis_state(inp, realAxisParameter)
    console["InitValues"] = [0 real(state.Z[1]) imag(state.Z[1]) real(state.delta[1]) imag(state.delta[1]) nothing]
    console = printTableHeader(console, log_file)

    β = 1 / (kb * itemp)
    Z_new = state.Z
    delta_new = state.delta
    data = [Z_new[1], delta_new[1]]

    # Iterate
    for i_it in 1:N_it
        delta_prev = copy(delta_new)
        Z_prev = copy(Z_new)
        broyden_beta = mixing_parameter(inp, i_it)
        gap0 = real(delta_prev[1])

        Z_new, delta_new = realEliashbergEq(muc_ME, β, delta_prev, Kp_func, Km_func, w_dynam, w_static, reOmega_c, i_it)

        Z_new = (1.0 - abs(broyden_beta)) .* Z_prev .+ abs(broyden_beta) .* Z_new
        delta_new = (1.0 - abs(broyden_beta)) .* delta_prev .+ abs(broyden_beta) .* delta_new

        convergence = sqrt(sum(abs2.(delta_new .- delta_prev))/length(delta_new))
        data = [Z_new[1], delta_new[1]]
        outputVec = real_axis_cdOS_output(i_it, Z_new, delta_new, convergence, gap0)
        print_real_axis_iteration(outputVec, console, log_file)

        if abs(convergence / gap0) < conv_thr && i_it > maximum([min_it, nItFullCoul + 1])
            print_real_axis_converged(itemp, console, log_file)

            state.Z = Z_new
            state.delta = delta_new
            if !vDOS_initial_guess
                save_real_axis_cDOS_outputs(itemp, inp, state, log_file)
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
            end

            return data, state
        end

        if i_it == N_it
            print_real_axis_not_converged(inp, console, log_file)
            data[2] = NaN

            # initial guess vDOS
            if  vDOS_initial_guess
                state.Z = Z_new
                state.delta = delta_new
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
    (_, _, dos_en, dos, Weep, dosef, idx_ef, _, _, _) = matval
    (; muc_ME, mu_flag, N_it, conv_thr, minGap, nItFullCoul, min_it, plot_flag) = inp
    (Kp_func, Km_func, w_static, _, w_static_chi) = realAxisParameter


    printTextCentered("T = " * string(itemp) * " K ", console["partingLine"], file=log_file, bold=true)
    printTee(log_file, "\n")

    # ----- electronic spectral information ----- #
    eps_1, eps_2 = transpose(dos_en[1:end-1]), transpose(dos_en[2:end])
    dos_1, dos_2 = transpose(dos[1:end-1]), transpose(dos[2:end])
    dε = eps_2 .- eps_1
    ddos = dos_2 .- dos_1
    electronic_spec = (eps_1, eps_2, dos_1, dos_2, dε, ddos)

    # ----- initial guesses -----#
    Z_new = state.Z
    delta_new = state.delta
    chi_new = state.chi
    phi_new = if inp.include_Weep == 1
            isnothing(state.phi) ? repeat(delta_new .* Z_new, 1, length(dos_en)) : state.phi
    else
        nothing
    end
    fermi_level = state.fermi_level

    # print console table
    console["InitValues"] = [0 real(state.Z[1]) imag(state.Z[1]) real(state.chi[1]) imag(state.chi[1]) state.fermi_level real(state.delta[1]) imag(state.delta[1]) nothing]
    console = printTableHeader(console, log_file)

    β = 1 / (kb * itemp)
    data = [state.Z[1], state.delta[1], state.chi[1]]

    s_flag = true

    # ------ integration grid and kernels ----- #
    # w_prime = make_vDOS_wprime_grid(max(0.1, real(state.delta[1])), inp)
    # Kernel_minus, Kernel_plus = evaluate_Kernels(w_static, w_prime, Km_func, Kp_func)
    # Kernel_plus_chi = evaluate_Kernels(w_static_chi, w_prime, Kp_func)
    # Kernels = (Kernel_minus, Kernel_plus, Kernel_plus_chi)

    for i_it in 1:N_it
        delta_prev = copy(delta_new)
        Z_prev = copy(Z_new)
        chi_prev = copy(chi_new)
        phi_prev = isnothing(phi_new) ? nothing : copy(phi_new)
        broyden_beta = mixing_parameter(inp, i_it)
        gap0 = real(delta_prev[1])
        w_prime = make_vDOS_wprime_grid(gap0, inp)
        wgCoulomb = minimum([1, i_it / nItFullCoul])


        if mu_flag == 1 
            @time fermi_level = mu_update_real_axis(itemp, fermi_level, w_static, w_static_chi, w_prime, dos_en, dos, Z_prev, delta_prev, chi_prev, inp.outdir, i_it)
        end

        if inp.include_Weep == 1
            @time Z_new, delta_new, chi_new, phi_new, s_flag = realEliashbergEq(β, Z_prev, phi_prev, chi_prev, Kp_func, Km_func, w_prime, w_static, w_static_chi, dosef, dos_en, dos, Weep, idx_ef, fermi_level, wgCoulomb, electronic_spec, i_it, gap0, s_flag, inp.outdir)
        else
            Z_new, delta_new, chi_new, s_flag = realEliashbergEq(muc_ME, β, Z_prev, delta_prev, chi_prev, Kp_func, Km_func, w_prime, w_static, w_static_chi, dosef, dos_en, dos, fermi_level, i_it, gap0, s_flag, inp.outdir)
        end

        chi_new = (1.0 - abs(broyden_beta)) .* chi_prev .+ abs(broyden_beta) .* chi_new
        Z_new = (1.0 - abs(broyden_beta)) .* Z_prev .+ abs(broyden_beta) .* Z_new
        delta_new = (1.0 - abs(broyden_beta)) .* delta_prev .+ abs(broyden_beta) .* delta_new
        if inp.include_Weep == 1
            phi_new = (1.0 - abs(broyden_beta)) .* phi_prev .+ abs(broyden_beta) .* phi_new
            delta_new = phi_new[:, idx_ef] ./ Z_new
        end

        convergence = sqrt(sum(abs2.(delta_new .- delta_prev))/length(delta_new))
        data = [Z_new[1], delta_new[1], chi_new[1]]
        outputVec = real_axis_vDOS_output(i_it, Z_new, delta_new, chi_new, fermi_level, convergence, gap0)
        print_real_axis_iteration(outputVec, console, log_file)

        if abs(convergence / gap0) < conv_thr && i_it > maximum([min_it, nItFullCoul + 1])
            print_real_axis_converged(itemp, console, log_file)

            if plot_flag
                selfEnergy = (delta_new, Z_new)
                plotSelfEnergyAtT(inp, itemp, selfEnergy, w_static)
                plotSelfEnergyAtT(inp, itemp, (chi_new,), w_static_chi, names=["chi"], labels=["χ(ω) / meV"])
            end

            # next guess
            state.Z = Z_new
            state.chi = chi_new
            state.delta = delta_new
            state.phi = phi_new
            state.fermi_level = fermi_level
            return data, state
        end

        if real(data[2]) < minGap && i_it > maximum([min_it, nItFullCoul + 1])
            print_real_axis_gap_too_small(itemp, minGap, console, log_file)
            data[2] = NaN
            return data, state
        end

        if i_it == N_it
            print_real_axis_not_converged(inp, console, log_file)
            data[2] = NaN
            return data, state
        end
    end
end


##############################################################
# -------------------- Helper functions -------------------- #
##############################################################
function initial_real_axis_state(inp, realAxisParameter)
    (; numReal_c) = inp
    (_, _, _, _, w_static_chi) = realAxisParameter

    return RealAxisState(
        ones(ComplexF64, numReal_c),
        ones(ComplexF64, numReal_c) .* (0.1 + im * 1e-4),
        -zeros(ComplexF64, length(w_static_chi)),
        0.0,
    )
end

"""
    mixing_parameter(inp, i_it)

linear mixing factor of eliashberg solutions
"""
function mixing_parameter(inp, i_it)
    if inp.mixing_beta == -1
        return maximum([0.5, 1.0 - 0.05 * (i_it - 1)])
    end

    return inp.mixing_beta
end


"""
    make_vDOS_wprime_grid(gap0, inp)

ω'-integration grid for vDOS
"""
function make_vDOS_wprime_grid(gap0, inp)
    (; num_wp1, num_wp2, wp_max, reOmega_c_shift,numReal_c_shift) = inp

    wp_inside = gap0 .* cos.((2 .* (0:num_wp1-1) .+ 1) ./ (2*num_wp1) .* π)
    wp_half = wp_max .+ wp_max .* cos.((2 .* (floor(num_wp2/2):num_wp2-1) .+ 1) ./ (2*num_wp2) .* π)
    wp_half = reverse(wp_half) .+ gap0
    wp_rest = range(maximum(wp_half), reOmega_c_shift, length=numReal_c_shift) 

    wp = vcat(
        reverse(wp_inside[wp_inside .> 0]),
        wp_half,
        collect(wp_rest),
    )

    return wp[wp .> 0]
end


function save_real_axis_cDOS_outputs(itemp, inp, state, log_file)
    plotSelfEnergyAtT(inp, itemp, (state.delta, state.Z))

    if inp.flag_writeSelfEnergy == 1
        try
            w_real = range(0, inp.reOmega_c, inp.numReal_c)
            saveSelfEnergyComponents(itemp, inp, w_real, state.delta, state.Z, mode="realAxis")
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
function initialize_real_axis_vDOS(itemp, inp, console, realAxisParameter, log_file)
    _, state = solve_realAxis_cDOS(itemp, inp, console, realAxisParameter, log_file; vDOS_initial_guess=true)
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
    outputVec, strConsole, format = formatTableRow(outputVec, console["width"], console["precision"])
    for i in axes(strConsole, 1)
        Printf.format(stdout, Printf.Format(strConsole[i]), format[i, 1], " ", format[i, 2], format[i, 3], outputVec[i], format[i, 4], " ")
    end

    for i in axes(strConsole, 1)
        Printf.format(log_file, Printf.Format(strConsole[i]), format[i, 1], " ", format[i, 2], format[i, 3], outputVec[i], format[i, 4], " ")
    end
end


"""
    print_real_axis_converged(itemp, console, log_file)

Print temperature converged
"""
function print_real_axis_converged(itemp, console, log_file)
    println(replace(console["Hline"], "." => " "))
    printstyled("\nConvergence achieved for T = " * string(itemp) * " K\n"; bold=false)

    println(log_file, replace(console["Hline"], "." => " "))
    printstyled(log_file, "\nConvergence achieved for T = " * string(itemp) * " K\n"; bold=false)
end


"""
    print_real_axis_gap_too_small(itemp, minGap, console, log_file)

gap at temperature too small
"""
function print_real_axis_gap_too_small(itemp, minGap, console, log_file)
    println(replace(console["Hline"], "." => " "))
    printstyled("\nTemperature (T = " * string(itemp) * " K) too high, gap value already smaller than " * string(round(minGap, digits=2)) * " meV!\n\n"; bold=false)

    println(log_file, replace(console["Hline"], "." => " "))
    printstyled(log_file, "\nTemperature (T = " * string(itemp) * " K) too high, gap value already smaller than " * string(round(minGap, digits=2)) * " meV!\n\n"; bold=false)
end


"""
     print_real_axis_not_converged(inp, console, log_file)

print max iterations exceeded
"""
function print_real_axis_not_converged(inp, console, log_file)
    println(replace(console["Hline"], "." => " "))
    printstyled("\nConvergence not achieved within " * string(inp.N_it) * " iterations\n"; bold=true)
    println("\n")

    println(log_file, replace(console["Hline"], "." => " "))
    printstyled(log_file, "\nConvergence not achieved within " * string(inp.N_it) * " iterations\n"; bold=true)
    println(log_file, "\n")
end

"""
    real_axis_cdOS_output(i_it, Z_new, delta_new, convergence, gap0)

cDOS console output
"""
function real_axis_cdOS_output(i_it, Z_new, delta_new, convergence, gap0)
    return [i_it, real(Z_new[1]), imag(Z_new[1]), real(delta_new[1]), imag(delta_new[1]), abs(convergence / gap0)]
end

"""
    real_axis_vDOS_output(i_it, Z_new, delta_new, chi_new, fermi_level, convergence, gap0)

vDOS console output
"""
function real_axis_vDOS_output(i_it, Z_new, delta_new, chi_new, fermi_level, convergence, gap0)
    return [i_it, real(Z_new[1]), imag(Z_new[1]), real(chi_new[1]), imag(chi_new[1]), fermi_level, real(delta_new[1]), imag(delta_new[1]), abs(convergence / gap0)]
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
