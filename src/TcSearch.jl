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
# inspired by the EPW implementation
# 2023-10-17 - Christoph Heil

export EliashbergSolver

"""
    solve_eliashberg(itemp, inp, console, matval)

Solve the eliashberg eq. self-consistently for a fixed temperature
"""
function solve_eliashberg(itemp, inp, console, matval, log_file)
    # destruct inputs
    (a2f_omega_fine, a2f_fine, dos_en, dos, Weep, dosef, idx_ef, ndos, BCS_gap, idxEncut) = matval
    (; cDOS_flag, include_Weep, omega_c, mixing_beta, nItFullCoul, muc_ME, mu_flag, outdir) = inp

    ### Matsubara frequencies ###
    beta = 1 / (kb * itemp)
    M = ceil(Int, (omega_c / (pi * kb * itemp) - 1) / 2)
    wsi = (2 * collect(1:M+1) .- 1) .* π .* kb .* itemp
    nsiw = size(wsi, 1)

    ### sparse sampling, consider only subset of mat frequencies up to M
    # only reasonable if T < 1 K
    if itemp < 5
        sparse_sampling_flag = 1
        ind_mat_freq = initSparseSampling(beta, omega_c, M)

        
        # sparse sampling, roman fbw paper
        #=
        fsparse = 1.0
        ind_mat_freq = cumsum(Int64.(round.(exp.(collect(1:nsiw)./(nsiw*fsparse)))))
        ind_mat_freq = ind_mat_freq[ind_mat_freq .<= nsiw]
        if ind_mat_freq[end] != nsiw
            ind_mat_freq = push!(ind_mat_freq, nsiw)
        end
        =#
        

        # write to console
        printTextCentered("T = "*string(itemp)*" K ", console["partingLine"], file = log_file, bold = true)
        printstyled("\n - Number of Matsubara Frequencies = ", length(ind_mat_freq), " / ", nsiw)
        println("\n")

        # write to log file
        printstyled(log_file, "\n - Number of Matsubara Frequencies = ", length(ind_mat_freq), " / ", nsiw)
        println(log_file, "\n")
    else
        sparse_sampling_flag = 0
        ind_mat_freq = collect(1:M+1)

        # write to console
        printTextCentered("T = "*string(itemp)*" K ", console["partingLine"], file = log_file, bold = true)
        printstyled("\n - Number of Matsubara Frequencies = ", nsiw)
        println("\n")

        # write to log file
        printstyled(log_file, "\n - Number of Matsubara Frequencies = ", nsiw)
        println(log_file, "\n")
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
            phici = -ones(ndos).*0.1
            phiphi = ones(nsiw) .* maximum([BCS_gap, 2*phici[1]])
            muintr = 0.0

            ### Print to console & log file
            console["InitValues"] = [0 phici[idx_ef] phiphi[1] znormi[1] shifti[1] -muintr deltai[idx_ef, 1] nothing]
            console = printTableHeader(console, log_file)

        elseif cDOS_flag == 1
            ### Initialize
            deltai = ones(ndos, nsiw) .* BCS_gap
            znormi = ones(nsiw) #.* 2 ./(1 .+ exp.(wsi .- round(nsiw/4))) .+1  # sigmoid function
            phici = -ones(ndos).*0.1
            phiphi = ones(nsiw) .* maximum([BCS_gap, 2*phici[1]])

            ### Print to console & log file
            console["InitValues"] = [0 phici[idx_ef] phiphi[1] znormi[1] deltai[idx_ef, 1] nothing]
            console = printTableHeader(console, log_file)

        end

    elseif include_Weep == 0

        if cDOS_flag == 0
            ### Initialize 
            deltai = ones(nsiw) .* BCS_gap
            znormi = ones(nsiw) #.* 2 ./(1 .+ exp.(wsi .- round(nsiw/4))) .+1  # sigmoid function
            shifti = zeros(nsiw)
            muintr = 0.0

            ### Print to console & log file
            console["InitValues"] = [0 znormi[1] shifti[1] -muintr deltai[1] nothing]
            console = printTableHeader(console, log_file)

        elseif cDOS_flag == 1
            ### Initialize 
            deltai = ones(nsiw) .* BCS_gap
            znormi = ones(nsiw) #.* 2 ./(1 .+ exp.(wsi .- round(nsiw/4))) .+1  # sigmoid function

            ### Print to console & log file
            console["InitValues"] = [0 znormi[1] deltai[1] nothing]
            console = printTableHeader(console, log_file)

        end

    else
        print(log_file, "ERROR: Unkwon mode! Check if the cDOS_flag and include_Weep flag are set correctly!")
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
                if mu_flag == 1 && i_it > 1
                    muintr = update_mu_own(itemp, wsi, dos_en, dos, znormip, deltaip, shiftip, idxEncut, outdir)
                end

                new_data = eliashberg_eqn(itemp, nsiw, wsi, ind_mat_freq, sparse_sampling_flag, lambdai, dosef, ndos, dos_en, dos, Weep, znormip, phiphip, phicip, shiftip, wgCoulomb, muintr, idxEncut)
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
                outputVec = [i_it, phici[idx_ef], phiphi[1], znormi[1], shifti[1], -muintr, deltai[idx_ef, 1], err_delta]

                # data for return
                data = [znormi[1], deltai[idx_ef, 1], shifti[1], -muintr]

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
                if mu_flag == 1 && i_it > 1
                    muintr = update_mu_own(itemp, wsi, dos_en, dos, znormip, deltaip, shiftip, idxEncut, outdir)
                end

                new_data = eliashberg_eqn(itemp, nsiw, wsi, ind_mat_freq, sparse_sampling_flag, lambdai, dos_en, dos, dosef, znormip, deltaip, shiftip, muc_ME, muintr, wgCoulomb, idxEncut)
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
                outputVec = [i_it, znormi[1], shifti[1], -muintr, deltai[1], err_delta]

                # data for return
                data = [znormi[1], deltai[1], shifti[1], -muintr]

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
        for i in axes(strConsole, 1)
           Printf.format(log_file, Printf.Format(strConsole[i]), format[i, 1], " ", format[i, 2], format[i, 3], outputVec[i], format[i, 4], " ")
        end



        ##### check convergence & termination criterion #####
        minIt = 15
        if err_delta < inp.conv_thr && i_it > maximum([minIt, inp.nItFullCoul+1])
            println(replace(console["Hline"], "." => " "))
            printstyled("\nConvergence achieved for T = " * string(itemp) * " K\n"; bold=false)

            println(log_file, replace(console["Hline"], "." => " "))
            printstyled(log_file, "\nConvergence achieved for T = " * string(itemp) * " K\n"; bold=false)

            # save Z, Delta, chi, phi
            if inp.flag_writeSelfEnergy == 1
                try 
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
                catch ex

                    # crash file
                    crashFile = open(inp.outdir * "CRASH", "a")
                    error_msg = sprint(showerror, ex)
                    st = sprint((io,v) -> show(io, "text/plain", v), stacktrace(catch_backtrace()))
                    print(crashFile, error_msg*"\n"*st)
                    print("\n\n")
                    close(crashFile)
        
                    # log file
                    print(log_file, "\n\nERROR while saving self energy components: ")
                    showerror(log_file, ex)
                    print(log_file, "\n\nFor further information please refer to the CRASH file\n\n")
        
                    # console
                    print(@red "\n\nERROR")
                    print(" while saving self energy components: ")
                    showerror(stdout, ex)
                    print("\n\nFor further information please refer to the CRASH file\n\n")
                end
            end
            

#=
            plot(wsi, shifti)
            savefig("shift_nomu_"*string(inp.Wcut)*".png")

            plot(wsi, znormip)
            savefig("Z_nomu_"*string(inp.Wcut)*".png")

            plot(wsi, deltai)
            savefig("delta_nomu_"*string(inp.Wcut)*".png")
=#
            return data
            break
        end

        # Gap too small
        if data[2] < 0.1 && i_it > maximum([minIt, inp.nItFullCoul+1])
            println(replace(console["Hline"], "." => " "))
            printstyled("\nTemperature (T = " * string(itemp) * " K) too high, gap value already smaller than 0.1 meV!\n\n"; bold=false)

            println(log_file, replace(console["Hline"], "." => " "))
            printstyled(log_file, "\nTemperature (T = " * string(itemp) * " K) too high, gap value already smaller than 0.1 meV!\n\n"; bold=false)

            data[2] = NaN
            return data
            break
        end

        # max number iterations reached
        if i_it == inp.N_it
            println(replace(console["Hline"], "." => " "))
            printstyled("\nConvergence not achieved within " * string(inp.N_it) * " iterations\n"; bold=true)
            println("\n")

            # log file
            println(log_file, replace(console["Hline"], "." => " "))
            printstyled(log_file, "\nConvergence not achieved within " * string(inp.N_it) * " iterations\n"; bold=true)
            println(log_file, "\n")
    

            data[2] = NaN
            return data
            break
        end


    end
end


# For each temperature solve eliashberg equations
function findTc(inp, console, matval, ML_Tc, log_file)
    inp.temps = sort(inp.temps)
    nT = size(inp.temps, 1)
    Delta0 = Vector{Float64}()
    Shift0 = Vector{Float64}()
    Znorm0 = Vector{Float64}()
    EfMu = Vector{Float64}()
    Tc = [NaN, NaN]


    if inp.TcSearchMode_flag == 0
        for iT in 1:nT
            itemp = inp.temps[iT]

            # solve Eliashberg equations
            data = solve_eliashberg(itemp, inp, console, matval, log_file)
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
            data = solve_eliashberg(itemp, inp, console, matval, log_file)
            if inp.cDOS_flag == 0
                Znorm0 = push!(Znorm0, data[1])
                Delta0 = push!(Delta0, data[2])
                Shift0 = push!(Shift0, data[3])
                EfMu = push!(EfMu, data[4])
            elseif inp.cDOS_flag == 1
                Znorm0 = push!(Znorm0, data[1])
                Delta0 = push!(Delta0, data[2])
            end


            # Escape
            if itemp < 1 && isnan(Delta0[end])
                # log file
                print(log_file, "\nTemperature already below 1 K. Material is most likely not a superconductor!\n")


                print("\nTemperature already below 1 K. Material is most likely not a superconductor!\n")

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
                    itemp = 1/2 
                end

            elseif sum(.~isnan.(Delta0)) == 1 
                # get a second gap value
                if length(Delta0) == 1
                    itemp = itemp + round(maximum([itemp / 3, 2]))
                else
                    itemp = ceil((maximum(inp.temps[.~isnan.(Delta0)]) + minimum(inp.temps[isnan.(Delta0)]))/2)
                end

            elseif sum(.~isnan.(Delta0)) == 2 && fitFlag
                # fit only once
                fitFlag = false

                nnanDelta = .~isnan.(Delta0)
                try
                    # fit gap values
                    p0 = convert(Vector{Float64}, [maximum(Delta0[nnanDelta]), 1, minimum(Delta0[nnanDelta])])
                    fit = curve_fit(m, inp.temps[nnanDelta], Delta0[nnanDelta], p0)
                    par = fit.param

                    # find root
                    m2(x) = m(x, par)
                    itemp = floor(find_zero(m2, par[2]))
                catch
                    itemp = maximum(inp.temps[nnanDelta]) + sum(inp.temps[nnanDelta]) / 2
                end

                # sanity check
                if length(Delta0) == 2  # no nans
                    if itemp < maximum(inp.temps)       # error in fit
                        itemp = ceil(maximum(inp.temps)*3/2)
                    elseif itemp == maximum(inp.temps)  # fit equals highest value
                        itemp += maximum([2, round(itemp/10)])
                    end
                elseif length(Delta0) >= 2
                    # prevent search below converged T, above not converged T
                    if itemp <= maximum(inp.temps[.~isnan.(Delta0)]) || itemp >= minimum(inp.temps[isnan.(Delta0)])
                        itemp = ceil((maximum(inp.temps[.~isnan.(Delta0)]) + minimum(inp.temps[isnan.(Delta0)]))/2)
                    end

                end
            else
                # search around fit value
                if any(isnan.(Delta0)) 
                    itemp = ceil((maximum(inp.temps[.~isnan.(Delta0)]) + minimum(inp.temps[isnan.(Delta0)]))/2)
                else
                    itemp += maximum([2, round(itemp/10)])
                end
            end
            
        end

    else
        error("Unknown Tc search mode! Please change the TcSearchMode_flag to an valid value!")
    end

    printTextCentered("Stopping now!", console["partingLine"], file = log_file, bold = true)

    return Tc, inp.temps, Znorm0, Delta0, Shift0, EfMu

end


"""
    EliashbergSolver(arguments)

Main function. User has to pass the input arguments and it returns the Tc.
"""
function EliashbergSolver(inp::arguments, testFlag=false)

    dt = @elapsed begin

        ### Create directory
        inp, log_file = createDirectory(inp::arguments)

        ### Check input
        try
            inp = checkInput(inp, log_file)
        catch ex
            crashFile = open(inp.outdir * "CRASH", "a")
            showerror(crashFile, ex)
            print("\n\n")
            close(crashFile)

            rethrow(ex)
        end

        ### read inputs
        inp, console, matval, ML_Tc = InputParser(inp, log_file)

        ### Print to console ###
        printFlagsAsText(inp, log_file)

        ########### start loop over temperatures ##########
        Tc = [NaN, NaN] 
        temps = Vector{Float64}()
        Delta0 = Vector{Float64}()
        Shift0 = Vector{Float64}()
        Znorm0 = Vector{Float64}()
        EfMu = Vector{Float64}()
        try
            Tc, temps, Znorm0, Delta0, Shift0, EfMu = findTc(inp, console, matval, ML_Tc, log_file)
        catch ex

            crashFile = open(inp.outdir * "CRASH", "a")
            print(crashFile, current_exceptions())
            print(crashFile, "\n\n")
            close(crashFile)

            print(log_file, "\nERROR: ")
            showerror(log_file, ex)
            print(log_file, "\n\nFor further information about the error please refer to the CRASH file\n")
            close(log_file)

            rethrow(ex)
        end
    

        ### write Tc to console
        try
            printSummary(inp, Tc, log_file)
        catch ex

            # crash file
            crashFile = open(inp.outdir * "CRASH", "a")
            print(crashFile, current_exceptions())
            print(crashFile, "\n\n")
            close(crashFile)

            # log file
            print(log_file, "\nERROR while printing the summary:\n")
            showerror(log_file, ex)
            print(log_file, "\n\nFor further information please refer to the CRASH file\n")

            # console
            print(@red "\nERROR")
            print(" while printing the summary:\n")
            showerror(stdout, ex)
            print("\nßnFor further information please refer to the CRASH file\n")
        end

        ### write output-file
        header = ["T / K", "Δ(0) / meV", "Z(0) / 1"]
        out_vars = zeros(size(Delta0, 1), 3)
        out_vars[:, 1] = temps
        out_vars[:, 2] = Delta0
        out_vars[:, 3] = Znorm0
        if inp.cDOS_flag == 0
            header = push!(header, "χ(0) / meV", "ϵ_F - μ / meV")
            out_vars = hcat(out_vars, Shift0, EfMu)
        end
        try
            writeToOutFile(Tc, inp, out_vars, header)
        catch ex

            # crash file
            crashFile = open(inp.outdir * "CRASH", "a")
            print(crashFile, current_exceptions())
            print(crashFile, "\n\n")
            close(crashFile)

            # log file
            print(log_file, "\nERROR while creating the summary file:\n")
            showerror(log_file, ex)
            print(log_file, "\n\nFor further information please refer to the CRASH file\n")

            # console
            print(@red "\nERROR")
            print(" while creating the summary file:\n\n")
            showerror(stdout, ex)
            #println(current_exceptions())  
            print("\n\nFor further information please refer to the CRASH file\n")

            #

        end

        ### figures
        if inp.flag_figure == 1
            try
                createFigures(inp, matval, Delta0, temps, log_file)
            catch ex
                
                # crash file
                crashFile = open(inp.outdir * "CRASH", "a")
                print(crashFile, current_exceptions())
                print(crashFile, "\n\n")
                close(crashFile)

                # log file
                print(log_file, "\nERROR while plotting: \n")
                showerror(log_file, ex)
                print(log_file, "\n\nFor further information please refer to the CRASH file\n")

                # console
                print(@red "\nERROR")
                print(" while plotting: ßn")
                showerror(stdout, ex)
                print("\n\nFor further information please refer to the CRASH file\n")
            end
        end
    end


    # print time elapsed
    print("\nTotal Runtime: ", dt, " seconds\n")
    # log file
    print(log_file, "\nTotal Runtime: ", dt, " seconds\n")

    # close & save
    close(log_file)


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





