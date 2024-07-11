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

export 
    findTc

function solve_eliashberg(itemp, gap0, console)
    """
    Solve eliashberg self-consistently

    -------------------------------------------------------------------
    Input:
        itemp:      Temperature

    --------------------------------------------------------------------
    Output:
        data:       Converged self energy on imaginary axis

    --------------------------------------------------------------------
    Comments:
        - 
    --------------------------------------------------------------------
    """


    ### Matsubara frequencies ###
    beta = 1 / (kb * itemp)
    M = ceil(Int, (omega_c / (pi * kb * itemp) - 1) / 2)
    wsi = (2 * collect(1:M+1) .- 1) .* π .* kb .* itemp
    nsiw = size(wsi, 1)

    ### sparse sampling, consider only subset of mat frequencies up to M
    # only reasonable if T < 1 K
    if itemp < 1
        global sparse_sampling_flag = 1
        ind_mat_freq = initSparseSampling(beta, omega_c, M)

        # write to console
        printTextCentered("T = "*string(itemp)*" K ", console["partingLine"] , true)
        printstyled("\n - Number of Matsubara Frequencies = ", length(ind_mat_freq), " / ", nsiw)
        println("\n")

        # write to log file
        if flag_log == 1
            printstyled(log_file, "\n - Number of Matsubara Frequencies = ", length(ind_mat_freq), " / ", nsiw)
            println(log_file, "\n")
        end
    else
        global sparse_sampling_flag = 0
        ind_mat_freq = collect(1:M+1)

        # write to console
        printTextCentered("T = "*string(itemp)*" K ", console["partingLine"] , true)
        printstyled("\n - Number of Matsubara Frequencies = ", nsiw)
        println("\n")

        # write to log file
        if flag_log == 1
            printstyled(log_file, "\n - Number of Matsubara Frequencies = ", nsiw)
            println(log_file, "\n")
        end
    end


    ### calc electron-phonon coupling ###
    lambdai = calcLambda(itemp, M, a2f_omega_fine, a2f_fine)


    ##### Initialize variables #####
    if include_Weep == 1
        if cDOS_flag == 0
            ### Initialize
            deltai = ones(ndos, nsiw) .* gap0
            znormi = ones(nsiw)
            shifti = -zeros(nsiw)
            phiphi = ones(nsiw) .* gap0
            phici = -zeros(ndos)
            muintr = ef

            ### Print to console & log file
            console["InitValues"] = [0 phici[idx_ef] phiphi[1] znormi[1] shifti[1] ef - muintr deltai[idx_ef, 1] nothing]
            console = printTableHeader(console)

        elseif cDOS_flag == 1
            ### Initialize
            deltai = ones(ndos, nsiw) .* gap0
            znormi = ones(nsiw)
            phiphi = ones(nsiw) .* gap0
            phici = zeros(ndos)

            ### Print to console & log file
            console["InitValues"] = [0 phici[idx_ef] phiphi[1] znormi[1] deltai[idx_ef, 1] nothing]
            console = printTableHeader(console)

        end

    elseif include_Weep == 0

        if cDOS_flag == 0
            ### Initialize 
            deltai = ones(nsiw) .* gap0
            znormi = ones(nsiw)
            shifti = zeros(nsiw)
            muintr = ef

            ### Print to console & log file
            console["InitValues"] = [0 znormi[1] shifti[1] ef - muintr deltai[1] nothing]
            console = printTableHeader(console)

        elseif cDOS_flag == 1
            ### Initialize 
            deltai = ones(nsiw) .* gap0
            znormi = ones(nsiw)

            ### Print to console & log file
            console["InitValues"] = [0 znormi[1] deltai[1] nothing]
            console = printTableHeader(console)

        end

    else
        @error "Unkwon mode! Check if the cDOS_flag and include_Weep flag are set correctly!"
    end




    ##### Start iterations #####
    err_delta = 0
    for i_it in 1:N_it
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
        broyden_beta = (1.0 - i_it / N_it)
        broyden_beta = maximum([0.5, 1.0 - 0.05*(i_it-1)])  

        if include_Weep == 1
            # weight weep
            wgWeep = minimum([1, i_it / nItFullWeep])

            ### Variable DoS ###
            if cDOS_flag == 0
                ### mu update 
                if mu_flag == 1
                    muintr = update_mu_own(itemp, wsi, M, muintr, dos_en, dos, znormip, deltaip, shiftip)
                end

                new_data = eliashberg_eqn(itemp, nsiw, wsi, ind_mat_freq, lambdai, dosef, ndos, dos_en, dos, znormip, phiphip, phicip, shiftip, wgWeep, muintr)
                shifti = (1.0 - abs(broyden_beta)) .* shifti .+ abs(broyden_beta) .* new_data[4]

                ### Constant DoS ###
            elseif cDOS_flag == 1
                new_data = eliashberg_eqn(itemp, nsiw, wsi, ind_mat_freq, lambdai, ndos, dos_en, dos, znormip, phiphip, phicip, wgWeep)
            end

            # linear mixing
            znormi = (1.0 - abs(broyden_beta)) .* znormi .+ abs(broyden_beta) .* new_data[1]
            phiphi = (1.0 - abs(broyden_beta)) .* phiphi .+ abs(broyden_beta) .* new_data[2]
            phici = (1.0 - abs(broyden_beta)) .* phici .+ abs(broyden_beta) .* new_data[3]
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
                    muintr = update_mu_own(itemp, wsi, M, muintr, dos_en, dos, znormip, deltaip, shiftip)
                end

                new_data = eliashberg_eqn(itemp, nsiw, wsi, ind_mat_freq, lambdai, dos_en, dos, dosef, znormip, deltaip, shiftip, muc_ME, muintr)
                shifti = (1.0 - abs(broyden_beta)) .* shifti .+ abs(broyden_beta) .* new_data[3]

            elseif cDOS_flag == 1
                new_data = eliashberg_eqn(itemp, nsiw, wsi, ind_mat_freq, lambdai, deltaip, muc_ME)
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
                Printf.format(log_file, Printf.Format(strConsole[i]), format[i, 1], " ", format[i, 2], format[i, 3], outputVec[i], format[i, 4], " ")
            end
        end


        ##### check convergence & termination criterion #####
        if err_delta < conv_thr
            println(replace(console["Hline"], "." => " "))
            printstyled("\nConvergence achieved for T = " * string(itemp) * " K\n"; bold=true)

            if flag_log == 1
                println(log_file, replace(console["Hline"], "." => " "))
                printstyled(log_file, "\nConvergence achieved for T = " * string(itemp) * " K\n\n"; bold=true)
            end

            return data
            break
        end
        # Gap too small
        if data[2] < 0.1 && i_it > 20
            println(replace(console["Hline"], "." => " "))
            printstyled("\nTemperature (T = " * string(itemp) * " K) too high, gap value already smaller than 0.1 meV!\n\n"; bold=true)

            data[2] = NaN
            return data
            break
        end
        # max number iterations reached
        if i_it == N_it
            println(replace(console["Hline"], "." => " "))
            printstyled("\nConvergence not achieved within " * string(N_it) * " iterations\n"; bold=true)
            println("\n")

            if flag_log == 1
                println(log_file, replace(console["Hline"], "." => " "))
                printstyled(log_file, "\nConvergence not achieved within " * string(N_it) * " iterations\n"; bold=true)
                println(log_file, "\n")
    
            end

            data[2] = NaN
            return data
            break
        end


    end
end


# Apparently loops in function are much faster than in global scope??
# Furthermore Julia does not like variables defined in global scope and soft scope
function findTc(temps, console)
    temps = sort(temps)
    nT = size(temps, 1)
    Delta0 = Vector{Float64}()
    Shift0 = Vector{Float64}()
    Znorm0 = Vector{Float64}()
    EfMu = Vector{Float64}()
    gap0 = 1#BCS_gap

    if TcSearchMode_flag == 0
        for iT in 1:nT
            itemp = temps[iT]

            # solve Eliashberg equations
            data = solve_eliashberg(itemp, gap0, console)
            if cDOS_flag == 0
                Znorm0 = push!(Znorm0, data[1])
                Delta0 = push!(Delta0, data[2])
                Shift0 = push!(Shift0, data[3])
                EfMu = push!(EfMu, data[4])
            elseif cDOS_flag == 1
                Znorm0 = push!(Znorm0, data[1])
                Delta0 = push!(Delta0, data[2])
            end

            if isnan(Delta0[end]) # this occurs if Delta0 < 0.1meV or max. number of iterations was exceeded
                # escape
                break
            end
        end

    elseif TcSearchMode_flag == 1

        # initial guess
        # rewrite s.t. user can specify array of temps which are all used for fitting
        if isempty(temps) || nT > 1
            itemp = round(AD_Tc)
        else
            if typeof(temps) == Array{Float64,1}
                itemp = temps[1]
            else
                itemp = temps
            end
        end

        # Init
        # expansion of a + b*log(c-x) at x = 0, use other function instead??
        m(x, p) = p[1] + log(p[3]) .- p[2] * x ./ p[3] .- p[2] * x .^ 2 / (2 * p[3]^2) .- p[2] * x .^ 3 / (3 * p[3]^3) .- p[2] * x .^ 4 / (4 * p[3]^4) .- p[2] * x .^ 5 / (5 * p[3]^5)
        temps = Vector{Float64}()
        fitFlag = true
        while true
            # save iterations
            temps = push!(temps, itemp)

            # solve Eliashberg equations
            data = solve_eliashberg(itemp, gap0, console)
            if cDOS_flag == 0
                Znorm0 = push!(Znorm0, data[1])
                Delta0 = push!(Delta0, data[2])
                Shift0 = push!(Shift0, data[3])
                EfMu = push!(EfMu, data[4])
            elseif cDOS_flag == 1
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

            elseif abs(temps[end] - temps[end-1]) <= 1 && ((isnan(Delta0[end]) && ~isnan(Delta0[end-1])) || (~isnan(Delta0[end]) && isnan(Delta0[end-1])))
                # converged
                break

            elseif sum(.~isnan.(Delta0)) == 2 && fitFlag
                # fit only once
                fitFlag = false

                # fit gap values
                nnanDelta = .~isnan.(Delta0)
                p0 = convert(Vector{Float64}, [maximum(Delta0[nnanDelta]), 1, minimum(Delta0[nnanDelta])])
                fit = curve_fit(m, temps[nnanDelta], Delta0[nnanDelta], p0)
                par = fit.param

                # find root
                m2(x) = m(x, par)
                itemp = floor(find_zero(m2, par[2]))

                # check if T < 0, happens when gap increases with T 
                # ideally a increasing gap should not happen at all
                if itemp < 0
                    itemp = maximum(temps[nnanDelta]) + sum(temps[nnanDelta]) / 2
                end

            else
                # search around fit value
                if isnan(Delta0[end]) 
                    itemp -= 1
                else
                    itemp += 1
                end
            end


            # Escape
            if itemp < 1
                if flag_log == 1
                    print(log_file, "\nTemperature already below 1 K. Material is most likely not a superconductor!\n")
                end

                print("\nTemperature already below 1 K. Material is most likely not a superconductor!\n")

                break
            elseif length(temps) > 500
                if flag_log == 1
                    print(log_file, "Couldn't find a Tc! \n")
                end

                print("Couldn't find a Tc! \n")

                break
            end
        end

    else
        error("Unknown Tc search mode! Please change the TcSearchMode_flag to an valid value!")
    end

    printTextCentered("Stopping now!", console["partingLine"], true)

    return temps, Znorm0, Delta0, Shift0, EfMu

end



############## main program starts here ##############
### defining constants ###
const Ry2meV = 13605.662285137
const THz2meV = 4.13566553853599;
const kb = 0.08617333262; # meV/K

### defining convergence parameters for Eliashberg solver
N_it = 5000;        # max. number of iterations for the Eliashberg solver
conv_thr = 1e-4;    # convergence threshold for Delta0 for solving the Eliashberg equations

### Check if input files are specified ###
cDOS_flag, include_Weep = checkInputFiles(cDOS_flag, include_Weep)

### Init table size ###
console = Dict()
if include_Weep == 1 && cDOS_flag == 0
    # header table
    console["header"] = ["it", "phic", "phiph", "znormi", "shifti", "ef-mu", "deltai", "err_delta"]
    # width table
    console["width"] = [8, 10, 10, 10, 10, 10, 10, 12]
    # precision data
    console["precision"] = [0, 2, 2, 2, 2, 2, 2, 5]

elseif include_Weep == 1 && cDOS_flag == 1
    # header table
    console["header"] = ["it", "phic", "phiph", "znormi", "deltai", "err_delta"]
    #width table
    console["width"] = [8, 10, 10, 10, 10, 12]
    # precision data
    console["precision"] = [0, 2, 2, 2, 2, 5]

elseif include_Weep == 0 && cDOS_flag == 0
    # header table
    console["header"] = ["it", "znormi", "shifti", "ef-mu", "deltai", "err_delta"]
    #width table
    console["width"] = [8, 10, 10, 10, 10, 12]
    # precision data
    console["precision"] = [0, 2, 2, 2, 2, 5]

elseif include_Weep == 0 && cDOS_flag == 1
    # header table
    console["header"] = ["it", "znormi", "deltai", "err_delta"]
    #width table
    console["width"] = [8, 14, 14, 14]
    # precision data
    console["precision"] = [0, 4, 4, 5]

else
    @error "Unkwon mode! Check if the cDOS_flag and include_Weep flag are set correctly!"
end
console = formatTableHeader(console)

console = printStartMessage(console)


########## READ-IN ##########
##### a2f file #####
a2f_omega_fine, a2f_fine = readIn_a2f(a2f_file, ind_smear, a2f_unit, nheader_a2f, nfooter_a2f, nsmear)

# determine superconducting properties from Allen-Dynes McMillan equation based on interpolated a2F
AD_data = calc_AD_Tc(a2f_omega_fine, a2f_fine, muc)
ML_Tc = AD_data[1] / kb;    # ML-Tc in K
AD_Tc = AD_data[2] / kb;    # AD-Tc in K
BCS_gap = AD_data[3];       # BCS gap value in meV
lambda = AD_data[4];        # total lambda
omega_log = AD_data[5];     # omega_log in meV

# convert muc for Migdal-Eliashberg
global muc_ME = muc / (1 + muc*log(200/omega_c))    
# CHANGE 200 to e.g. maximum(a2f_omega_fine[a2f_fine .> 0.1])

# print Allen-Dynes
printADtable(console)


########## read-in and process Dos and Weep ##########
### Weep ###
if include_Weep == 1
    Weep, unitWeepFile = readIn_Weep(Weep_file, Weep_unit, nheader_Weep, nfooter_Weep)

    # Include full weep at iteration:
    nItFullWeep = 5
    # if greater than 10 adapt termination criterion for min iterations !!
end


### Dos ###
if @isdefined(dosW_file) && @isdefined(dos_file) && include_Weep == 1
    # QE/Abinit/... dos
    dos_en, dos, ef, unitDosFile = readIn_Dos(dos_file, colFermi_dos, spinDos, dos_unit, nheader_dos, nfooter_dos)
    # W dos
    dosW_en, dosW, ef, unitDosWFile = readIn_Dos(dosW_file, colFermi_dosW, spinDosW, dosW_unit, nheader_dosW, nfooter_dosW)

    # overlap of energies
    en_interval = [dosW_en[findfirst(dosW_en .> dos_en[1])]; dosW_en[findlast(dosW_en .< dos_en[end])]]
    # number of points overlapping
    idxOverlap = [findfirst(dosW_en .> dos_en[1]), findlast(dosW_en .< dos_en[end])]
    Nitp = idxOverlap[2] - idxOverlap[1] + 1
    
    # restrict Weep
    Weep = Weep[idxOverlap[1]:idxOverlap[2], idxOverlap[1]:idxOverlap[2]]
    
    # Interpolation Dos
    dos_en, dos = interpolateDos(dos_en, dos, en_interval, Nitp)

elseif @isdefined(dos_file)
     # QE/Abinit/... dos
     dos_en, dos, ef, unitDosFile = readIn_Dos(dos_file, colFermi_dos, spinDos, dos_unit, nheader_dos, nfooter_dos)

elseif @isdefined(dosW_file) 
     # W dos
     dos_en, dos, ef, unitDosFile = readIn_Dos(dosW_file, colFermi_dosW, spinDosW, dosW_unit, nheader_dosW, nfooter_dosW)
 
end


### REMOVE - START ###
# restrict Dos/Weep to a subset of grid points
# ONLY RELEVANT TO TEST THE SCALING OF THE CODE
if ~isnothing(Nrestrict)
    if include_Weep == 1
        # call restrict function
        Weep, dos, dos_en = restrictInput(Nrestrict, wndRestrict, ef, dos, dos_en, Weep)
    else
        # call restrict function
        dos, dos_en = restrictInput(Nrestrict, wndRestrict, ef, dos, dos_en)
    end
end
### REMOVE - END ###


### remove zeros at begining/end of dos ###
if include_Weep == 1 
    dos, dos_en, Weep = neglectZeros(dos, dos_en, Weep)
else
    dos, dos_en = neglectZeros(dos, dos_en)
end


### length energy vector ###
ndos = size(dos_en, 1)

### index of fermi energy ###
idx_ef = findmin(abs.(dos_en .- ef)) # find value closest to ef
idx_ef = idx_ef[2] # index of closest value

### dos at ef ###
dosef = dos[idx_ef]


###### Print to console #####
# Flags
printFlagsAsText()


########### start loop over temperatures ##########
temps, Znorm0, Delta0, Shift0, EfMu = findTc(temps, console)

# write Tc to console
printSummary()

# write output-file
if flag_outfile == 1
    out_var = zeros(nT, 5)
    out_var[:, 1] = temps
    out_var[:, 2] = Delta0
    if cDOS_flag == 0
        out_var[:, 3] = Znorm0
        out_var[:, 4] = Shift0
        out_var[:, 5] = EfMu
    elseif cDOS_flag == 1
        out_var[:, 3] = Znorm0
    end

    if mu_flag == 0
        writedlm(outdir * "/Delta0_constmu_sm" * string(ind_smear) * ".txt", out_var)
    elseif mu_flag == 1
        writedlm(outdir * "/Delta0_varmu_sm" * string(ind_smear) * ".txt", out_var)
    end
end

# print a2F vs. energy
if flag_figure == 1
    plot_font = "Computer Modern"
    default(
        fontfamily=plot_font,
        linewidth=2,
        framestyle=:box,
        label=nothing,
        grid=false
    )

    xlim_max = round(maximum(a2f_omega_fine) / 10 * 1.01, RoundUp) * 10
    xtick_val = 0:10:xlim_max
    ylim_max = round(maximum(a2f_fine), RoundUp)

    plot(a2f_omega_fine, a2f_fine)
    xlims!(0, xlim_max)
    ylims!(0, ylim_max)
    title!(material)
    xlabel!(L"\omega ~ \mathrm{(meV)}")
    ylabel!(L"\alpha^2F ~ \mathrm{(1/meV)}")
    savefig(outdir * "/a2F_sm" * string(ind_smear) * ".pdf")
end

# print gap vs. temperature
if flag_figure == 1
    nan_ind = Delta0 .== Delta0
    Delta0 = Delta0[nan_ind]
    temps = temps[nan_ind]

    plot_font = "Computer Modern"
    default(
        fontfamily=plot_font,
        linewidth=2,
        framestyle=:box,
        label=nothing,
        grid=false
    )

    if maximum(temps) < 10
        xlim_max = round(maximum(temps) * 1.1, RoundUp)
        xtick_val = 0:1:xlim_max
    else
        xlim_max = round(maximum(temps) / 10 * 1.01, RoundUp) * 10
        xtick_val = 0:10:xlim_max
    end

    if maximum(Delta0) < 10
        ylim_max = round(maximum(Delta0), RoundUp)
    else
        ylim_max = round(maximum(Delta0) / 10, RoundUp) * 10
    end

    plot(temps, Delta0, marker=:circle)
    plot!(xticks=xtick_val)
    xlims!(0, xlim_max)
    ylims!(0, ylim_max)
    title!(material)
    xlabel!(L"T ~ \mathrm{(K)}")
    ylabel!(L"\Delta_0 ~ \mathrm{(meV)}")

    if cDOS_flag == 0
        if mu_flag == 0
            savefig(outdir * "/Delta0_constmu_mu" * string(muc) * "_sm" * string(ind_smear) * ".pdf")
        elseif mu_flag == 1
            savefig(outdir * "/Delta0_varmu_mu" * string(muc) * "_sm" * string(ind_smear) * ".pdf")
        end
    else
        savefig(outdir * "/Delta0_conDOS_mu" * string(muc) * "_sm" * string(ind_smear) * ".pdf")
    end
end