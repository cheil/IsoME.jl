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


"""
    solve_eliashberg(itemp, inp, console, matval, log_file)

Solve the eliashberg eq. self-consistently for a fixed temperature
"""
function solve_eliashberg(itemp, inp, console, matval, log_file)
    # destruct inputs
    (a2f_omega_fine, a2f_fine, dos_en, dos, Weep, dosef, idx_ef, ndos, BCS_gap, idxShiftcut) = matval
    (; cDOS_flag, include_Weep, imOmega_c, mixing_beta, nItFullCoul, muc_ME, mu_flag, outdir, sparseSamplingTemp, flag_acon) = inp

    # active output table (mode 0 populates cDOS or vDOS by cDOS_flag)
    table = cDOS_flag == 1 ? console.cDOS : console.vDOS

    ### Matsubara frequencies ###
    beta = 1 / (kb * itemp)
    M = ceil(Int, (imOmega_c / (pi * kb * itemp) - 1) / 2)
    wsi = collect((2 * (0:M) .+ 1) .* π .* kb .* itemp)
    nsiw = size(wsi, 1)

    ### sparse sampling
    if itemp < sparseSamplingTemp    
        sparse_sampling_flag = 1
        ind_mat_freq = initSparseSampling(beta, imOmega_c, M)       

        # write to console
        printTextCentered("T = "*string(itemp)*" K ", console.partingLine, file = log_file, bold = true)
        printstyled("\n - Number of Matsubara Frequencies = ", length(ind_mat_freq), " / ", nsiw)
        println("\n")

        # write to log file
        printstyled(log_file, "\n - Number of Matsubara Frequencies = ", length(ind_mat_freq), " / ", nsiw)
        println(log_file, "\n")
    else
        sparse_sampling_flag = 0
        ind_mat_freq = collect(1:M+1)

        # write to console
        printTextCentered("T = "*string(itemp)*" K ", console.partingLine, file = log_file, bold = true)
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
            znormi = ones(nsiw) 
            shifti = -zeros(nsiw)
            phici = -ones(ndos).*0.1
            phiphi = ones(nsiw) .* maximum([BCS_gap, 2*phici[1]])
            fermi_level = 0.0

            ### Print to console & log file
            initValues = [0 phici[idx_ef] phiphi[1] znormi[1] shifti[1] -fermi_level deltai[idx_ef, 1] nothing]
            printTableHeader(table, initValues, log_file)

        elseif cDOS_flag == 1
            ### Initialize
            deltai = ones(ndos, nsiw) .* BCS_gap
            znormi = ones(nsiw) 
            phici = -ones(ndos).*0.1
            phiphi = ones(nsiw) .* maximum([BCS_gap, 2*phici[1]])

            ### Print to console & log file
            initValues = [0 phici[idx_ef] phiphi[1] znormi[1] deltai[idx_ef, 1] nothing]
            printTableHeader(table, initValues, log_file)

        end

    elseif include_Weep == 0

        if cDOS_flag == 0
            ### Initialize 
            deltai = ones(nsiw) .* BCS_gap
            znormi = ones(nsiw) 
            shifti = zeros(nsiw)
            fermi_level = 0.0

            ### Print to console & log file
            initValues = [0 znormi[1] shifti[1] -fermi_level deltai[1] nothing]
            printTableHeader(table, initValues, log_file)

        elseif cDOS_flag == 1
            ### Initialize 
            deltai = ones(nsiw) .* BCS_gap
            znormi = ones(nsiw) 

            ### Print to console & log file
            initValues = [0 znormi[1] deltai[1] nothing]
            printTableHeader(table, initValues, log_file)

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
        if isnan(mixing_beta)
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
                    fermi_level = update_mu_own(itemp, wsi, dos_en, dos, znormip, deltaip, shiftip, idxShiftcut, fermi_level, outdir)
                end

                new_data = eliashberg_eqn(itemp, nsiw, wsi, ind_mat_freq, sparse_sampling_flag, lambdai, dosef, ndos, dos_en, dos, Weep, znormip, phiphip, phicip, shiftip, wgCoulomb, fermi_level, idxShiftcut)
                shifti = (1.0 - abs(broyden_beta)) .* shiftip .+ abs(broyden_beta) .* new_data[4]

                ### Constant DoS ###
            elseif cDOS_flag == 1
                new_data = eliashberg_eqn(itemp, nsiw, wsi, ind_mat_freq, sparse_sampling_flag, lambdai, ndos, dos_en, dos, Weep, znormip, phiphip, phicip, idx_ef, wgCoulomb)
            end

            # linear mixing
            znormi = (1.0 - abs(broyden_beta)) .* znormip .+ abs(broyden_beta) .* new_data[1]
            phiphi = (1.0 - abs(broyden_beta)) .* phiphip .+ abs(broyden_beta) .* new_data[2]
            phici  = (1.0 - abs(broyden_beta)) .* phicip  .+ abs(broyden_beta) .* new_data[3]
            deltai = (phiphi' .+ phici) ./ znormi'

            rel_delta = sum(abs.(deltai[idx_ef, :] .- deltaip[idx_ef, :]))
            abs_delta = sum(abs.(deltai[idx_ef, :]))
            err_delta = rel_delta / abs_delta


            ### Console Output ###
            if cDOS_flag == 0
                # Console output
                outputVec = [i_it, phici[idx_ef], phiphi[1], znormi[1], shifti[1], -fermi_level, deltai[idx_ef, 1], err_delta]

                # data for return
                data = [znormi[1], deltai[idx_ef, 1], shifti[1], -fermi_level]

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
                    fermi_level = update_mu_own(itemp, wsi, dos_en, dos, znormip, deltaip, shiftip, idxShiftcut, fermi_level, outdir)
                end

                new_data = eliashberg_eqn(itemp, nsiw, wsi, ind_mat_freq, sparse_sampling_flag, lambdai, dos_en, dos, dosef, znormip, deltaip, shiftip, muc_ME, fermi_level, wgCoulomb, idxShiftcut)
                shifti = (1.0 - abs(broyden_beta)) .* shiftip .+ abs(broyden_beta) .* new_data[3]

            elseif cDOS_flag == 1
                new_data = eliashberg_eqn(itemp, nsiw, wsi, ind_mat_freq, sparse_sampling_flag, lambdai, deltaip, muc_ME, wgCoulomb)
            end

            znormi = (1.0 - abs(broyden_beta)) .* znormip .+ abs(broyden_beta) .* new_data[1]
            deltai = (1.0 - abs(broyden_beta)) .* deltaip .+ abs(broyden_beta) .* new_data[2]

            rel_delta = sum(abs.(deltai - deltaip))
            abs_delta = sum(abs.(deltai))
            err_delta = rel_delta / abs_delta


            ### Console Output ###
            if cDOS_flag == 0
                # Console output
                outputVec = [i_it, znormi[1], shifti[1], -fermi_level, deltai[1], err_delta]

                # data for return
                data = [znormi[1], deltai[1], shifti[1], -fermi_level]

            elseif cDOS_flag == 1
                # Console output
                outputVec = [i_it, znormi[1], deltai[1], err_delta]

                # data for return
                data = [znormi[1], deltai[1]]

            end # cDOS_flag

        end # include_Weep


        ##### Print to console #####
        outputVec, strConsole, format = formatTableRow(outputVec, table.width, table.precision)
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
            println(replace(table.Hline, "." => " "))
            printstyled("\nConvergence achieved for T = " * string(itemp) * " K\n"; bold=false)

            println(log_file, replace(table.Hline, "." => " "))
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
                    writeToCrashFile(inp)

                    # console / log file
                    printWarning("Error while saving self energy components.", log_file, ex=ex)
                end
            end

            ### Analytic Continuation ###
            if flag_acon
                try
                    if include_Weep == 1
                        if inp.cDOS_flag == 0
                            acon(inp, itemp, wsi, nsiw, deltai, znormi, log_file, idx_ef=idx_ef, shifti = shifti)
                        elseif inp.cDOS_flag == 1
                            acon(inp, itemp, wsi, nsiw, deltai, znormi, log_file, idx_ef=idx_ef)
                        end
                    elseif inp.include_Weep == 0
                        if inp.cDOS_flag == 0
                            acon(inp, itemp, wsi, nsiw, deltai, znormi, log_file, shifti = shifti)                     
                        elseif inp.cDOS_flag == 1
                            acon(inp, itemp, wsi, nsiw, deltai, znormi, log_file)          
                        end
                    end
                catch ex
                    # crash file
                    writeToCrashFile(inp)

                    # console / log file
                    printWarning("Error in analytic continuation.", log_file, ex=ex)
                end
            end 


            return data
            break
        end

        # Gap too small
        if data[2] < inp.minGap && i_it > maximum([minIt, inp.nItFullCoul+1])
            println(replace(table.Hline, "." => " "))
            printstyled("\nTemperature (T = " * string(itemp) * " K) too high, gap value already smaller than "*string(round(inp.minGap, digits=2))*" meV!\n\n"; bold=false)

            println(log_file, replace(table.Hline, "." => " "))
            printstyled(log_file, "\nTemperature (T = " * string(itemp) * " K) too high, gap value already smaller than "*string(round(inp.minGap, digits=2))*" meV!\n\n"; bold=false)

            data[2] = NaN
            return data
            break
        end

        # max number iterations reached
        if i_it == inp.N_it
            println(replace(table.Hline, "." => " "))
            printstyled("\nConvergence not achieved within " * string(inp.N_it) * " iterations\n"; bold=true)
            println("\n")

            # log file
            println(log_file, replace(table.Hline, "." => " "))
            printstyled(log_file, "\nConvergence not achieved within " * string(inp.N_it) * " iterations\n"; bold=true)
            println(log_file, "\n")
    

            data[2] = NaN
            return data
            break
        end


    end
end

 
"""
    findTc(inp, console, matval, ML_Tc, log_file)

Start the Tc search mode or solve the imaginary eliashberg equations for each temperature
"""
function findTc(inp, console, matval, ML_Tc, log_file)
    inp.temps = sort(inp.temps)
    nT = size(inp.temps, 1)
    Delta0 = Vector{Float64}()
    Shift0 = Vector{Float64}()
    Znorm0 = Vector{Float64}()
    EfMu = Vector{Float64}()
    Tc = [NaN, NaN]


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

            # solve Eliashberg equations
            data = solve_eliashberg(itemp, inp, console, matval, log_file)
            if inp.cDOS_flag == 0
                Znorm0 = push!(Znorm0, data[1])
                Delta0 = push!(Delta0, data[2])
                Shift0 = push!(Shift0, data[3])
                EfMu   = push!(EfMu, data[4])
            elseif inp.cDOS_flag == 1
                Znorm0 = push!(Znorm0, data[1])
                Delta0 = push!(Delta0, data[2])
            end


            # Escape
            if itemp < 1 && isnan(Delta0[end])
                # log file
                print(log_file, "Lowest temperature of Tc search mode reached. If you want to search at even lower temperatures consider setting them manually!\n")

                print("Lowest temperature of Tc search mode reached. If you want to search at even lower temperatures consider setting them manually!\n")

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
                    itemp = 1/2 
                end

            elseif sum(.~isnan.(Delta0)) == 1 
                # get a second gap value
                if length(Delta0) == 1
                    itemp = itemp + round(maximum([itemp / 2, 2]))      
                else
                    itemp = ceil((maximum(inp.temps[.~isnan.(Delta0)]) + minimum(inp.temps[isnan.(Delta0)]))/2)
                end

            elseif sum(.~isnan.(Delta0)) == 2 && fitFlag
                # fit only once
                fitFlag = false

                nnanDelta = .~isnan.(Delta0)
                # fit gap values
                p0 = convert(Vector{Float64}, [maximum(Delta0[nnanDelta]), 1, minimum(Delta0[nnanDelta])])
                try
                    fit = curve_fit(m, inp.temps[nnanDelta], Delta0[nnanDelta], p0)
                    par = fit.param

                    # find root
                    m2(x) = m(x, par)
                    itemp = floor(find_zero(m2, par[2]))

                catch ex
                    ex isa InterruptException && rethrow(ex)
                    # expansion to third order --> analytical formula for root (only one real root)
                    a=p0[1]
                    c=p0[3]
                    itemp = -c / 2
                    itemp += -(3 * c^2) / (2 * (5 * c^3 + 12 * a * c^3 + 
                                2 * sqrt(13 * c^6 + 30 * a * c^6 + 36 * a^2 * c^6))^(1/3))
                    itemp += 0.5 * (5 * c^3 + 12 * a * c^3 + 
                                2 * sqrt(13 * c^6 + 30 * a * c^6 + 36 * a^2 * c^6))^(1/3)
                    itemp = round(itemp)
                    # formula works only if a,c are far away from the Tc
                    # if a,c close to Tc the estimated T will be too small but this case is caputred by the sanity check
                end

                # sanity check
                if length(Delta0) == 2  # no nans
                    if itemp < maximum(inp.temps)       # error in fit
                        itemp = ceil(maximum(inp.temps)*3/2)
                    elseif itemp == maximum(inp.temps)  # fit equals highest converged value
                        itemp += maximum([2, round(itemp/10)])
                    end
                elseif length(Delta0) >= 2
                    # prevent search below converged T, above not converged T
                    if itemp <= maximum(inp.temps[nnanDelta]) || itemp >= minimum(inp.temps[isnan.(Delta0)])
                        itemp = ceil((maximum(inp.temps[nnanDelta]) + minimum(inp.temps[isnan.(Delta0)]))/2)
                    end

                end
            else
                # search around fit value
                if any(isnan.(Delta0)) 
                    itemp = ceil((maximum(inp.temps[.~isnan.(Delta0)]) + minimum(inp.temps[isnan.(Delta0)]))/2)
                else
                    itemp += maximum([2, round(itemp/5)])
                end
            end
            
        end

    else    # given temperatures
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

    end

    printTextCentered("Stopping now!", console.partingLine, file = log_file, bold = true)

    return Tc, inp.temps, Znorm0, Delta0, Shift0, EfMu

end


"""
    EliashbergSolver(inp)

Main function. User has to pass the input arguments and it returns the Tc.
"""
function EliashbergSolver(inp::arguments)

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
                checkInput(inp)
            end

            ### read inputs
            inp, console, matval, ML_Tc = stage("while reading the inputs") do
                InputParser(inp, log_file)
            end

            ### Print to console ###
            printFlagsAsText(inp, log_file)

            ########### start loop over temperatures ##########
            Tc, temps, Znorm0, Delta0, Shift0, EfMu = stage("while solving the Eliashberg equations") do
                findTc(inp, console, matval, ML_Tc, log_file)
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
                    header = "# T/K   Δ(0)/meV   Z(0)/1   "
                    out_vars = zeros(size(Delta0, 1), 3)
                    out_vars[:, 1] = temps
                    out_vars[:, 2] = Delta0
                    out_vars[:, 3] = Znorm0
                    if inp.cDOS_flag == 0
                        header = header * "χ(0)/meV   ϵ_F-μ/meV   "
                        out_vars = hcat(out_vars, Shift0, EfMu)
                    end
                    createSummaryFile(inp, Tc, out_vars, header)
                end


                ### figures
                if inp.flag_figure == 1
                    attempt(inp, log_file, "Error while plotting. Skipping plots.") do
                        createFigures(inp, matval, Delta0, temps, Tc, log_file)
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





