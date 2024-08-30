#
#
#  _                 __  __   ______ 
# (_)               |  \/  | |  ____|
#  _   ___    ___   | \  / | | |__   
# | | / __|  / _ \  | |\/| | |  __|  
# | | \__ \ | (_) | | |  | | | |____ 
# |_| |___/  \___/  |_|  |_| |______|
#                                    
#                                    
#
#
# routine to solve the full-bandwidth isotropic Migdal-Eliashberg equations
# adapted from the EPW implementation
# 2023-10-17 - Christoph Heil

export EliashbergSolver

"""
    solve_eliashberg(itemp, inp, console, matval)

Solve the eliashberg eq. self-consistently for a fixed temperature
"""
function solve_eliashberg(itemp, inp, console, matval)
    # destruct inputs
    (a2f_omega_fine, a2f_fine, dos_en, dos, Weep, ef, dosef, idx_ef, ndos, BCS_gap) = matval
    (; cDOS_flag, include_Weep, flag_log, omega_c, mixing_beta, nItFullCoul, muc_ME, mu_flag) = inp

    ### Matsubara frequencies ###
    beta = 1 / (kb * itemp)
    M = ceil(Int, (omega_c / (pi * kb * itemp) - 1) / 2)
    wsi = (2 * collect(1:M+1) .- 1) .* π .* kb .* itemp
    nsiw = size(wsi, 1)

    ### sparse sampling, consider only subset of mat frequencies up to M
    # only reasonable if T < 1 K
    if itemp < 5
        sparse_sampling_flag = 1
        #ind_mat_freq = initSparseSampling(beta, omega_c, M)

        # sparse sampling, roman fbw paper
        fsparse = 1.0
        ind_mat_freq = cumsum(Int64.(round.(exp.(collect(1:nsiw)./(nsiw*fsparse)))))
        ind_mat_freq = ind_mat_freq[ind_mat_freq .<= nsiw]
        if ind_mat_freq[end] != nsiw
            ind_mat_freq = push!(ind_mat_freq, nsiw)
        end

        # write to console
        printTextCentered(inp, "T = "*string(itemp)*" K ", console["partingLine"] , true)
        printstyled("\n - Number of Matsubara Frequencies = ", length(ind_mat_freq), " / ", nsiw)
        println("\n")

        # write to log file
        if flag_log == 1
            printstyled(inp.log_file, "\n - Number of Matsubara Frequencies = ", length(ind_mat_freq), " / ", nsiw)
            println(inp.log_file, "\n")
        end
    else
        sparse_sampling_flag = 0
        ind_mat_freq = collect(1:M+1)

        # write to console
        printTextCentered(inp, "T = "*string(itemp)*" K ", console["partingLine"] , true)
        printstyled("\n - Number of Matsubara Frequencies = ", nsiw)
        println("\n")

        # write to log file
        if flag_log == 1
            printstyled(inp.log_file, "\n - Number of Matsubara Frequencies = ", nsiw)
            println(inp.log_file, "\n")
        end
    end


    ### calc electron-phonon coupling ###
    lambdai = calcLambda(itemp, M, a2f_omega_fine, a2f_fine)


    ##### Initialize variables #####
    if include_Weep == 1
        if cDOS_flag == 0
            ### Initialize
            deltai = ones(ndos, nsiw) .* BCS_gap
            znormi = ones(nsiw) #.* 2 ./(1 .+ exp.(wsi .- round(nsiw/4))) .+1  # sigmoid function
            shifti = -zeros(nsiw)
            phiphi = ones(nsiw) .* BCS_gap
            phici = -zeros(ndos)
            muintr = ef

            ### Print to console & log file
            console["InitValues"] = [0 phici[idx_ef] phiphi[1] znormi[1] shifti[1] ef - muintr deltai[idx_ef, 1] nothing]
            console = printTableHeader(inp, console)

        elseif cDOS_flag == 1
            ### Initialize
            deltai = ones(ndos, nsiw) .* BCS_gap
            znormi = ones(nsiw) #.* 2 ./(1 .+ exp.(wsi .- round(nsiw/4))) .+1  # sigmoid function
            phiphi = ones(nsiw) .* BCS_gap
            phici = zeros(ndos)

            ### Print to console & log file
            console["InitValues"] = [0 phici[idx_ef] phiphi[1] znormi[1] deltai[idx_ef, 1] nothing]
            console = printTableHeader(inp, console)

        end

    elseif include_Weep == 0

        if cDOS_flag == 0
            ### Initialize 
            deltai = ones(nsiw) .* BCS_gap
            znormi = ones(nsiw) #.* 2 ./(1 .+ exp.(wsi .- round(nsiw/4))) .+1  # sigmoid function
            shifti = zeros(nsiw)
            muintr = ef

            ### Print to console & log file
            console["InitValues"] = [0 znormi[1] shifti[1] ef - muintr deltai[1] nothing]
            console = printTableHeader(inp, console)

        elseif cDOS_flag == 1
            ### Initialize 
            deltai = ones(nsiw) .* BCS_gap
            znormi = ones(nsiw) #.* 2 ./(1 .+ exp.(wsi .- round(nsiw/4))) .+1  # sigmoid function

            ### Print to console & log file
            console["InitValues"] = [0 znormi[1] deltai[1] nothing]
            console = printTableHeader(inp, console)

        end

    else
        @error "Unkwon mode! Check if the cDOS_flag and include_Weep flag are set correctly!"
    end




    ##### Start iterations #####
    err_delta = 0
    for i_it in 1:inp.N_it
        if include_Weep == 1
            if cDOS_flag == 0
                deltaip = copy(deltai)
                znormip = znormi
                shiftip = shifti
                phiphip = phiphi
                phicip = phici
            elseif cDOS_flag == 1
                deltaip = copy(deltai)
                znormip = znormi
                phiphip = phiphi
                phicip = phici
            end
        elseif include_Weep == 0
            if cDOS_flag == 0
                deltaip = copy(deltai)
                znormip = znormi
                shiftip = shifti
            elseif cDOS_flag == 1
                deltaip = copy(deltai)
                znormip = znormi
            end
        end


        # mixing beta
        if mixing_beta == -1
            broyden_beta = maximum([0.5, 1.0 - 0.05*(i_it-1)]) 
        else
            broyden_beta = mixing_beta
        end 

        # weight coulomb interaction
        wgCoulomb = minimum([1, i_it / nItFullCoul])

        if include_Weep == 1
            ### Variable DoS ###
            if cDOS_flag == 0
                ### mu update 
                if mu_flag == 1
                    muintr = update_mu_own(itemp, wsi, M, ef, dos_en, dos, znormip, deltaip, shiftip)
                end

                new_data = eliashberg_eqn(itemp, nsiw, wsi, ind_mat_freq, sparse_sampling_flag, lambdai, dosef, ndos, dos_en, dos, Weep, znormip, phiphip, phicip, shiftip, wgCoulomb, muintr)
                shifti = (1.0 - abs(broyden_beta)) .* shifti .+ abs(broyden_beta) .* new_data[4]

                ### Constant DoS ###
            elseif cDOS_flag == 1
                new_data = eliashberg_eqn(itemp, nsiw, wsi, ind_mat_freq, sparse_sampling_flag, lambdai, ndos, dos_en, dos, Weep, znormip, phiphip, phicip, idx_ef, wgCoulomb)
            end

            # linear mixing
            znormi = (1.0 - abs(broyden_beta)) .* znormi .+ abs(broyden_beta) .* new_data[1]
            phiphi = (1.0 - abs(broyden_beta)) .* phiphi .+ abs(broyden_beta) .* new_data[2]
            phici  = (1.0 - abs(broyden_beta)) .* phici  .+ abs(broyden_beta) .* new_data[3]
            deltai = (phiphi' .+ phici) ./ znormi'

            rel_delta = sum(abs.(deltai[idx_ef, :] .- deltaip[idx_ef, :]))
            abs_delta = sum(abs.(deltai[idx_ef, :]))
            err_delta = rel_delta / abs_delta


            ### Console Output ###
            if cDOS_flag == 0
                # Console output
                outputVec = [i_it, phici[idx_ef], phiphi[1], znormi[1], shifti[1], ef - muintr, deltai[idx_ef, 1], err_delta]

                # data for return
                data = [znormi[1], deltai[idx_ef, 1], shifti[1], ef - muintr]

            elseif cDOS_flag == 1
                # Console output
                outputVec = [i_it, phici[idx_ef], phiphi[1], znormi[1], deltai[idx_ef, 1], err_delta]

                # data for return
                data = [znormi[1], deltai[idx_ef, 1]]

            end # cDOS_flag


        ##### No Weep #####
        elseif include_Weep == 0

            if cDOS_flag == 0
                ### mu update
                if mu_flag == 1
                    muintr = update_mu_own(itemp, wsi, M, ef, dos_en, dos, znormip, deltaip, shiftip)
                end

                new_data = eliashberg_eqn(itemp, nsiw, wsi, ind_mat_freq, sparse_sampling_flag, lambdai, dos_en, dos, dosef, znormip, deltaip, shiftip, muc_ME, muintr, wgCoulomb)
                shifti = (1.0 - abs(broyden_beta)) .* shifti .+ abs(broyden_beta) .* new_data[3]

            elseif cDOS_flag == 1
                new_data = eliashberg_eqn(itemp, nsiw, wsi, ind_mat_freq, sparse_sampling_flag, lambdai, deltaip, muc_ME, wgCoulomb)
            end

            znormi = (1.0 - abs(broyden_beta)) .* znormi .+ abs(broyden_beta) .* new_data[1]
            deltai = (1.0 - abs(broyden_beta)) .* deltai .+ abs(broyden_beta) .* new_data[2]

            rel_delta = sum(abs.(deltai - deltaip))
            abs_delta = sum(abs.(deltai))
            err_delta = rel_delta / abs_delta


            ### Console Output ###
            if cDOS_flag == 0
                # Console output
                outputVec = [i_it, znormi[1], shifti[1], ef - muintr, deltai[1], err_delta]

                # data for return
                data = [znormi[1], deltai[1], shifti[1], ef - muintr]

            elseif cDOS_flag == 1
                # Console output
                outputVec = [i_it, znormi[1], deltai[1], err_delta]

                # data for return
                data = [znormi[1], deltai[1]]

            end # cDOS_flag

        end # include_Weep


        ##### Print to console #####
        outputVec, strConsole, format = formatTableRow(outputVec, console["width"], console["precision"])
        for i in axes(strConsole, 1)
            Printf.format(stdout, Printf.Format(strConsole[i]), format[i, 1], " ", format[i, 2], format[i, 3], outputVec[i], format[i, 4], " ")
        end

        ### print to log file ###
        if flag_log == 1
            for i in axes(strConsole, 1)
                Printf.format(inp.log_file, Printf.Format(strConsole[i]), format[i, 1], " ", format[i, 2], format[i, 3], outputVec[i], format[i, 4], " ")
            end
        end


        ##### check convergence & termination criterion #####
        if err_delta < inp.conv_thr && i_it > maximum([15, inp.nItFullCoul+1])
            println(replace(console["Hline"], "." => " "))
            printstyled("\nConvergence achieved for T = " * string(itemp) * " K\n"; bold=false)

            if flag_log == 1
                println(inp.log_file, replace(console["Hline"], "." => " "))
                printstyled(inp.log_file, "\nConvergence achieved for T = " * string(itemp) * " K\n\n"; bold=false)
            end

            # save Z, Delta, chi, phi
            if inp.flag_writeSelfEnergy == 1
                if include_Weep == 1
                    if inp.cDOS_flag == 0
                        saveSelfEnergyComponents(itemp, inp, wsi, deltai, znormi, epsilon=dos_en, chi=shifti, phiph=phiphi, phic=phici)
                    elseif inp.cDOS_flag == 1
                        saveSelfEnergyComponents(itemp, inp, wsi, deltai, znormi, epsilon=dos_en, phiph=phiphi, phic=phici)
                    end
                elseif inp.include_Weep == 0
                    if inp.cDOS_flag == 0
                        saveSelfEnergyComponents(itemp, inp, wsi, deltai, znormi, chi=shifti)                           
                    elseif inp.cDOS_flag == 1
                        saveSelfEnergyComponents(itemp, inp, wsi, deltai, znormi)                
                    end
                end
            end

            return data
            break
        end

        # Gap too small
        if data[2] < 0.1 && i_it > maximum([15, inp.nItFullCoul+1])
            println(replace(console["Hline"], "." => " "))
            printstyled("\nTemperature (T = " * string(itemp) * " K) too high, gap value already smaller than 0.1 meV!\n\n"; bold=false)

            if flag_log == 1
                println(inp.log_file, replace(console["Hline"], "." => " "))
                printstyled(inp.log_file, "\nTemperature (T = " * string(itemp) * " K) too high, gap value already smaller than 0.1 meV!\n\n"; bold=false)
            end

            data[2] = NaN
            return data
            break
        end

        # max number iterations reached
        if i_it == inp.N_it
            println(replace(console["Hline"], "." => " "))
            printstyled("\nConvergence not achieved within " * string(inp.N_it) * " iterations\n"; bold=true)
            println("\n")

            if flag_log == 1
                println(inp.log_file, replace(console["Hline"], "." => " "))
                printstyled(inp.log_file, "\nConvergence not achieved within " * string(inp.N_it) * " iterations\n"; bold=true)
                println(inp.log_file, "\n")
    
            end

            data[2] = NaN
            return data
            break
        end


    end
end


# For each temperature solve eliashberg equations
function findTc(inp, console, matval, ML_Tc)
    inp.temps = sort(inp.temps)
    nT = size(inp.temps, 1)
    Delta0 = Vector{Float64}()
    Shift0 = Vector{Float64}()
    Znorm0 = Vector{Float64}()
    EfMu = Vector{Float64}()

    if inp.TcSearchMode_flag == 0
        for iT in 1:nT
            itemp = inp.temps[iT]

            # solve Eliashberg equations
            data = solve_eliashberg(itemp, inp, console, matval)
            if inp.cDOS_flag == 0
                Znorm0 = push!(Znorm0, data[1])
                Delta0 = push!(Delta0, data[2])
                Shift0 = push!(Shift0, data[3])
                EfMu = push!(EfMu, data[4])
            elseif inp.cDOS_flag == 1
                Znorm0 = push!(Znorm0, data[1])
                Delta0 = push!(Delta0, data[2])
            end

            if isnan(Delta0[end]) 
                # escape
                inp.temps = inp.temps[1:length(Delta0)]
                break
            end
        end

    elseif inp.TcSearchMode_flag == 1

        # initial guess, Machine learning Tc           
        itemp = round(ML_Tc)
        # rewrite s.t. user can specify array of temps which are all used for fitting ?

        # expansion of a + b*log(c-x) at x = 0, use other function instead??
        m(x, p) = p[1] + log(p[3]) .- p[2] * x ./ p[3] .- p[2] * x .^ 2 / (2 * p[3]^2) .- p[2] * x .^ 3 / (3 * p[3]^3) .- p[2] * x .^ 4 / (4 * p[3]^4) .- p[2] * x .^ 5 / (5 * p[3]^5)
        inp.temps = Vector{Float64}()
        fitFlag = true
        while true
            # save iterations
            inp.temps = push!(inp.temps, itemp)

            # solve Eliashberg equations
            data = solve_eliashberg(itemp, inp, console, matval)
            if inp.cDOS_flag == 0
                Znorm0 = push!(Znorm0, data[1])
                Delta0 = push!(Delta0, data[2])
                Shift0 = push!(Shift0, data[3])
                EfMu = push!(EfMu, data[4])
            elseif inp.cDOS_flag == 1
                Znorm0 = push!(Znorm0, data[1])
                Delta0 = push!(Delta0, data[2])
            end


            if all(isnan.(Delta0))
                # temperature too high
                itemp = floor(itemp / 2)

            elseif sum(.~isnan.(Delta0)) == 1 #|| sum(.~isnan.(Delta0)) == 2
                # get a second gap value
                if isnan(Delta0[end])
                    itemp = floor(itemp / 2)
                else
                    itemp = floor(itemp + maximum([itemp / 2, 2]))
                end

            elseif abs(inp.temps[end] - inp.temps[end-1]) <= 1 && ((isnan(Delta0[end]) && ~isnan(Delta0[end-1])) || (~isnan(Delta0[end]) && isnan(Delta0[end-1])))
                # converged
                break

            elseif sum(.~isnan.(Delta0)) == 2 && fitFlag
                # fit only once
                fitFlag = false

                # fit gap values
                nnanDelta = .~isnan.(Delta0)
                p0 = convert(Vector{Float64}, [maximum(Delta0[nnanDelta]), 1, minimum(Delta0[nnanDelta])])
                fit = curve_fit(m, inp.temps[nnanDelta], Delta0[nnanDelta], p0)
                par = fit.param

                # find root
                m2(x) = m(x, par)
                itemp = floor(find_zero(m2, par[2]))

                # check if T < 0, happens when gap increases with T 
                # ideally a increasing gap should not happen at all
                if itemp < 0
                    itemp = maximum(inp.temps[nnanDelta]) + sum(inp.temps[nnanDelta]) / 2
                end

            else
                # search around fit value
                if isnan(Delta0[end]) 
                    itemp = ceil((maximum(inp.temps[.~isnan.(Delta0)]) + inp.temps[end])/2)
                else
                    itemp += 1
                end
            end


            # Escape
            if itemp < 1
                if inp.flag_log == 1
                    print(inp.log_file, "\nTemperature already below 1 K. Material is most likely not a superconductor!\n")
                end

                print("\nTemperature already below 1 K. Material is most likely not a superconductor!\n")

                break
            elseif length(inp.temps) > 500
                if inp.flag_log == 1
                    print(inp.log_file, "Couldn't find a Tc! \n")
                end

                print("Couldn't find a Tc! \n")

                break
            end
        end

    else
        error("Unknown Tc search mode! Please change the TcSearchMode_flag to an valid value!")
    end

    printTextCentered(inp, "Stopping now!", console["partingLine"], true)

    return inp.temps, Znorm0, Delta0, Shift0, EfMu

end

"""
    EliashbergSolver(arguments)

Main function. User has to pass the input arguments and it returns the Tc.
"""
function EliashbergSolver(inp::arguments, testFlag = false)
    dt = @elapsed begin
        # inputs
        inp, console, matval, ML_Tc = InputParser(inp)
        
        ### Print to console ###
        printFlagsAsText(inp)

        ########### start loop over temperatures ##########
        temps, Znorm0, Delta0, Shift0, EfMu = findTc(inp, console, matval, ML_Tc)

        # write Tc to console
        printSummary(inp, temps, Delta0)

        # write output-file
        if inp.flag_outfile == 1
            ### Data ###
            header = ["T / K", "Δ(0) / meV", "Z(0) / 1"]

            out_vars = zeros(size(Delta0, 1), 3)
            out_vars[:, 1] = temps
            out_vars[:, 2] = Delta0
            out_vars[:, 3] = Znorm0
            if inp.cDOS_flag == 0
                header = push!(header, "χ(0) / meV", "ϵ_F - μ / meV")
                out_vars = hcat(out_vars, Shift0, EfMu)
            end

            writeToOutFile(inp, out_vars, header)

        end


        ### figures
        a2f_omega_fine, a2f_fine = matval
        plot_font = "Computer Modern"
        default(
            fontfamily=plot_font,
            linewidth=2,
            framestyle=:box,
            label=nothing,
            grid=false
        )
        if inp.flag_figure == 1
            # print a2F vs. energy
            xlim_max = round(maximum(a2f_omega_fine) / 10 * 1.01, RoundUp) * 10
            xtick_val = 0:10:xlim_max
            ylim_max = round(maximum(a2f_fine), RoundUp)

            plot(a2f_omega_fine, a2f_fine)
            xlims!(0, xlim_max)
            ylims!(0, ylim_max)
            title!(inp.material)
            xlabel!(L"\omega ~ \mathrm{(meV)}")
            ylabel!(L"\alpha^2F ~ \mathrm{(1/meV)}")
            savefig(inp.outdir * "/a2F_sm" * string(inp.ind_smear) * ".pdf")

            if all(isnan.(Delta0))
                print(@blue "Info: ")
                println("No superconducting gap found - skipping plot\n")
            else
                # print gap vs. temperature
                Delta0_plot = Delta0
                Delta0_plot[isnan.(Delta0)] .= 0
                temps_plot = temps

                plot_font = "Computer Modern"

                if maximum(temps_plot) < 10
                    xlim_max = round(maximum(temps_plot) * 1.1, RoundUp)
                    xtick_val = 0:1:xlim_max
                else
                    xlim_max = round(maximum(temps_plot) / 10 * 1.01, RoundUp) * 10
                    xtick_val = 0:10:xlim_max
                end

                if maximum(Delta0_plot) < 10
                    ylim_max = round(maximum(Delta0_plot), RoundUp)
                else
                    ylim_max = round(maximum(Delta0_plot) / 10, RoundUp) * 10
                end

                plot(temps_plot, Delta0_plot, marker=:circle)
                plot!(xticks=xtick_val)
                xlims!(0, xlim_max)
                ylims!(0, ylim_max)
                title!(inp.material)
                xlabel!(L"T ~ \mathrm{(K)}")
                ylabel!(L"\Delta_0 ~ \mathrm{(meV)}")

                namePlot = "Delta0"
                if inp.material != "material"
                    namePlot = namePlot*"_"*inp.material
                end
                if inp.cDOS_flag == 1
                    namePlot = namePlot*"_"*"cDOS"
                else
                    namePlot = namePlot*"_"*"vDOS"
                end
                if inp.include_Weep== 1
                    namePlot = namePlot*"_"*"W"
                else
                    namePlot = namePlot*"_"*"muc"
                end
                namePlot = namePlot*".pdf"
                savefig(inp.outdir*namePlot)
            end
        end
    end

    # print time elapsed
    print("Total Runtime: ", dt, " seconds\n")
    if inp.flag_log == 1
        print(inp.log_file, "Total Runtime: ", dt, " seconds\n")
   
        # close & save
        close(inp.log_file)
    end

    ### !!! Better solution for runtest !!!
    if testFlag
        if all(isnan.(Delta0))
            Tc = NaN
        else
            Tc = maximum(temps[.~isnan.(Delta0)])
        end
        return Tc
    end
end





