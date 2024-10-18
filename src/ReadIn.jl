"""
    File containing everything needed for the read in of the input files
        - alpha2F
        - Dos
        - Weep

    Julia Packages:
        - DelimitedFiles

    Comments:
        - 

"""


"""
    InputParser()

Include the file specified via path.

# Examples
"""
function InputParser(inp::arguments, log_file::IOStream)

    ### Init table size ###
    console = Dict()
    if inp.include_Weep == 1 && inp.cDOS_flag == 0
        # header table
        console["header"] = ["it", "phic", "phiph", "znormi", "shifti", "ef-mu", "deltai", "err_delta"]
        # width table
        console["width"] = [8, 10, 10, 10, 10, 10, 10, 12]
        # precision data
        console["precision"] = [0, 2, 2, 2, 2, 2, 2, 5]

    elseif inp.include_Weep == 1 && inp.cDOS_flag == 1
        # header table
        console["header"] = ["it", "phic", "phiph", "znormi", "deltai", "err_delta"]
        #width table
        console["width"] = [8, 10, 10, 10, 10, 12]
        # precision data
        console["precision"] = [0, 2, 2, 2, 2, 5]

    elseif inp.include_Weep == 0 && inp.cDOS_flag == 0
        # header table
        console["header"] = ["it", "znormi", "shifti", "ef-mu", "deltai", "err_delta"]
        #width table
        console["width"] = [8, 10, 10, 11, 10, 12]
        # precision data
        console["precision"] = [0, 2, 2, 2, 2, 5]

    elseif inp.include_Weep == 0 && inp.cDOS_flag == 1
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

    console = printStartMessage(console, log_file)
  

    ########## READ-IN ##########
    ####### a2f file #######
    a2f_omega, a2f = readIn_a2f(inp.a2f_file, inp.ind_smear, inp.a2f_unit, inp.nheader_a2f, inp.nfooter_a2f, inp.nsmear)


    ########## Dos and Weep ##########
    if ~isempty(inp.dos_file) 
        # read dos
        dos_en, dos, ef, inp.dos_unit = readIn_Dos(inp.dos_file, inp.ef, inp.colFermi_dos, inp.spinDos, inp.dos_unit, inp.nheader_dos, inp.nfooter_dos)

        # remove zeros at begining/end of dos
        #dos, dos_en = discardZeros(dos, dos_en)

        # Wcut in meV
        if inp.Wcut != -1
            if (inp.Wcut) > dos_en[end] || (inp.Wcut) < dos_en[1]
                print(@yellow "WARNING: ")
                println("Wcut larger than given energy interval!")
            else
                logWcut = (dos_en .> (inp.Wcut)) .& (dos_en .< (inp.Wcut))
                dos      = dos[logWcut]
                dos_en   = dos_en[logWcut]
            end
        end


        if inp.include_Weep == 1
            # read Weep + energy grid points
            Weep, inp.Weep_unit = readIn_Weep(inp.Weep_file, inp.Weep_col, inp.Weep_unit, inp.nheader_Weep, inp.nfooter_Weep)
            W_en = readIn_Wen(inp.Wen_file, inp.efW, inp.colFermi_Wen, inp.Wen_unit, inp.nheader_Wen, inp.nfooter_Wen)  

            ### Interpolate Weep & Dos onto non-uniform grid
            bndItp = [1000, 500]
            en_range = minimum([abs(dos_en[1]), abs(dos_en[end])])
            en_interval = [W_en[findfirst(W_en .> -en_range)], -bndItp[1], -bndItp[2], bndItp[2], bndItp[1],  W_en[findlast(W_en .< en_range)]]  
            gridSpecs = [("step", 50), ("step", 5), ("step", 1), ("step", 5), ("step", 50)]
                        
            # check if energy window < bndItp
            logBnd = (en_interval[1] .<= en_interval) .& (en_interval[end] .>= en_interval)
            en_interval = en_interval[logBnd]
            gridSpecs   = gridSpecs[append!(logBnd[2:3], [true], logBnd[4:5])]
          
            # interpolate 
            dos_en, dos = interpolateDos(dos_en, dos, en_interval, gridSpecs)
            Weep = interpolateWeep(W_en, Weep, en_interval, gridSpecs)

            # idx encut
            idxEncut = [findfirst(dos_en .> -inp.encut), findlast(dos_en .< inp.encut)]

        else
            # no Weep
            Weep = nothing
            W_en = nothing

            # Interpolation Dos
            bndItp = [1000, 500]
            en_range = minimum([abs(dos_en[1]), abs(dos_en[end])])
            en_interval = [-en_range, -bndItp[1], -bndItp[2], bndItp[2], bndItp[1],  en_range]  
            gridSpecs = [("step", 50), ("step", 5), ("step", 1), ("step", 5), ("step", 50)] 
            # check if energy window < bndItp
            logBnd = (en_interval[1] .<= en_interval) .& (en_interval[end] .>= en_interval)
            en_interval = en_interval[logBnd]
            gridSpecs   = gridSpecs[append!(logBnd[2:3], [true], logBnd[4:5])]
          
            # interpolate 
            dos_en, dos = interpolateDos(dos_en, dos, en_interval, gridSpecs)

            # idx encut
            idxEncut = [findfirst(dos_en .> -inp.encut), findlast(dos_en .< inp.encut)]
        end
    else
        dos = []
        dos_en = []
        ef = nothing
        idxEncut = -1
    end

    ### calc mu*'s
    if inp.include_Weep == 0
        if inp.mu != -1
            if isnothing(ef) || ef == -1 || ef == 0
                print(@yellow "WARNING: ")
                print("Unable to calculate μ* from μ without the fermi-energy! Either provide a dos-file or set the fermi-energy directly in the input structure.\n\n")
            
                # log file
                print(log_file, "WARNING: Unable to calculate μ* from μ without the fermi-energy! Either provide a dos-file or set the fermi-energy directly in the input structure.\n\n")
            else
                # ensure μ*_ME < 4*μ --> reasonable ??
                if inp.omega_c > ef*exp(3/(4*inp.mu))
                    inp.omega_c = ef*exp(3/(4*inp.mu))

                    print(@yellow "WARNING: ")
                    print("Matsubara cutoff would lead to μ*_ME > 4*μ. omega_c has been set to a smaller value.\n\n")
            
                    # log file
                    print(log_file, "WARNING: Matsubara cutoff would lead to μ*_ME > 4*μ. omega_c has been set to a smaller value.\n\n")
                end

                inp.muc_AD = inp.mu / (1 + inp.mu*log(ef/maximum(a2f_omega[a2f .> 0.01])))
                inp.muc_ME = inp.mu / (1 + inp.mu*log(ef/inp.omega_c))
            end
        elseif inp.muc_ME == -1
            inp.muc_ME = inp.muc_AD / (1 + inp.muc_AD*log(maximum(a2f_omega[a2f .> 0.01])/inp.omega_c))

            if inp.muc_ME < 0 || inp.muc_ME > 0.8 || inp.muc_ME > 3*inp.muc_AD                 
                inp.muc_ME = minimum([3*inp.muc_AD, 0.8])

                print(@yellow "WARNING: ")
                print("Couldn't calculate a reasonable μ*_ME from μ*_AD.\n Using μ*_ME = minimum(3*μ*_AD, 0.8) instead.\n Consider setting it manually! \n\n")
            
                # log file
                print(log_file, "WARNING: Couldn't find a reasonable μ*_ME from μ*_AD.\n Using μ*_ME = minimum(3*μ*_AD, 0.8) instead.\n Consider setting it manually! \n\n")
            end
        end
    end

    # restrict Dos/Weep to a subset of grid points
    if inp.Nrestrict != -1
        if inp.include_Weep == 1
            # call restrict function
            Weep, dos, dos_en = restrictInput(inp.Nrestrict, inp.wndRestrict, dos, dos_en, Weep)
        else
            # call restrict function
            dos, dos_en = restrictInput(inp.Nrestrict, inp.wndRestrict, dos, dos_en)
        end
    end


    if inp.include_Weep == 0 && inp.cDOS_flag == 1
        # default values in case no dos-file is given
        ndos = -1
        idx_ef = -1
        dosef = -1
    else
        # length energy vector
        ndos = size(dos_en, 1)

        # index of fermi energy
        idx_ef = findmin(abs.(dos_en))
        idx_ef = idx_ef[2]

        # dos at ef
        dosef = dos[idx_ef]
    end

    ### determine superconducting properties from Allen-Dynes McMillan equation based on interpolated a2F
    AD_data = calc_AD_Tc(a2f_omega, a2f, inp.muc_AD)
    ML_Tc = AD_data[1] / kb    # ML-Tc in K
    AD_Tc = AD_data[2] / kb    # AD-Tc in K
    BCS_gap = AD_data[3]       # BCS gap value in meV
    lambda = AD_data[4]        # total lambda
    omega_log = AD_data[5]     # omega_log in meV

    # print Allen-Dynes
    printADtable(console, ML_Tc, AD_Tc, BCS_gap, lambda, omega_log, log_file)
    

    # material specific values
    matval = (a2f_omega, a2f, dos_en, dos, Weep, dosef, idx_ef, ndos, BCS_gap, idxEncut)

    return inp, console, matval, ML_Tc
end


"""
    createDirectory()

Check which input files (a2f, dos, weep) exist. Overwrite flags if neccessary.
"""
function createDirectory(inp::arguments)

    # create output directory
    if isempty(inp.outdir) 
        inp.outdir = "./"
    elseif ~(inp.outdir[end] == '/' || inp.outdir[end] == '\\')
        inp.outdir = inp.outdir*"/"
    end

    idxDir = 1
    tempDir = inp.outdir[1:end-1]
    while isdir(inp.outdir)
        inp.outdir = tempDir*"_"*string(idxDir)*"/"
        idxDir+=1
    end

    try
        mkdir(inp.outdir)
    catch ex
        error("Couldn't write into " * inp.outdir * "! Outdir may not writable or an invalid path.\n\n"*"Exception: "*ex.msg)
    end

    # open log file
    if isfile(inp.outdir * "log.txt")
        rm(inp.outdir * "log.txt")
    end

    log_file = open(inp.outdir * "log.txt", "w")
    print(log_file)

    return inp, log_file

end


"""
    checkInput!()

Check which input files (a2f, dos, weep) exist. Overwrite flags if neccessary.
"""
function checkInput(inp::arguments, log_file::IOStream)

    # check input files / cDOS & Weep
    if ~isfile(inp.a2f_file)
        error("Invalid path to a2f-file!")
        
    elseif ~isfile(inp.dos_file) && (inp.cDOS_flag == 0 || inp.include_Weep == 1)
        inp.cDOS_flag = 1 
        inp.include_Weep = 0
        print(@yellow "WARNING: ")
        print("No valid Dos-file specified! Calculating Tc within constant Dos approximation using Anderson-pseudopotential instead\n\n")
    
        # log file
        print(log_file, "WARNING: No valid Dos-file specified! Calculating Tc within constant Dos approximation using Anderson-pseudopotential instead\n\n")
    
    elseif (~(isfile(inp.Weep_file) || ~isfile(inp.Wen_file))) && inp.include_Weep == 1
        inp.include_Weep = 0
        print(@yellow "WARNING: ")
        print("No valid Weep-files specified! Calculating Tc using Anderson-pseudopotential instead\n\n")
    
        # log file
        print(log_file, "WARNING: No valid Weep-files specified! Calculating Tc using Anderson-pseudopotential instead\n\n")
    end


    # Tc search mode
    if inp.temps == [-1] && inp.TcSearchMode_flag == 0
        inp.TcSearchMode_flag = 1
        print(@yellow "WARNING: ")
        print("No temperatures specified - Activating Tc search mode instead!\n\n")
    
        # log file
        print(log_file, "WARNING: No temperatures specified - Activating Tc search mode instead!\n\n")
    end

    return inp

end

"""
    readIn_a2f(a2f_file, indSmear [,unit, nheader, nfooter, nsmear])   

Read in a2f file to solve the isotropic Migdal-Eliashberg equations

The first column must contain the energies, the second column onwards a2F values for different smearings
"""
function readIn_a2f(a2f_file, indSmear=-1, unit="", nheader=-1, nfooter=-1, nsmear=-1)   
    ### Read in a2f file ###
    a2f_data = readdlm(a2f_file);

    ### Define defaults
    (nheader != -1) || (nheader = findfirst(isa.(a2f_data[:,1], Number))-1)
    (nfooter != -1) || (nfooter = size(a2f_data, 1) - findlast(isa.(a2f_data[:,1], Number)))
    (nsmear != -1) || (nsmear = length(a2f_data[nheader+1, isa.(a2f_data[nheader+1,:], Number)])-1)
    (indSmear != -1) || (indSmear = Int64(ceil(nsmear/2)))

    ### Remove header & footer
    header = join(a2f_data[1:nheader,:], " ")
    a2f_data = Float64.(a2f_data[nheader+1:end-nfooter, 1:nsmear+1])  # previous version: nheader+1:end-nfooter

    ### Convert omega ###
    omega_raw = a2f_data[:, 1]
    (~isempty(unit)) || (unit = getUnit(header, "a2F"))
    if "meV" == unit
        omega_raw = omega_raw
    elseif "eV" == unit
        omega_raw = omega_raw*1000
    elseif "THz" == unit
            omega_raw = omega_raw * THz2meV
    elseif "Ry" == unit
        omega_raw = omega_raw * Ry2meV
    else
        error("Invalid Unit! Please check the header of the a2F-file and try again!")
    end

    ### a2f for one smearing ###
    a2f_raw = a2f_data[:, indSmear+1] 

    ### interpolate a2F on 10x finer grid ###
    omega_fine = range(1e-2, stop=omega_raw[end], length=10 * size(omega_raw)[1])
    a2f_int = linear_interpolation(omega_raw, a2f_raw, extrapolation_bc=Line())
    a2f_fine = a2f_int(omega_fine)
    a2f_fine[a2f_fine.<0.0] .= 0.0


    return omega_fine, a2f_fine

end


# Read in DOS
function readIn_Dos(dos_file, ef =-1, colFermi=1, spin=2, unit="", nheader=-1, nfooter=-1, outdir ="")
    """
    Read in DoS file to solve the isotropic Migdal-Eliashberg equations
    All quantities are converted to meV

    -------------------------------------------------------------------
    Input:
        Dos_file:   path to the Dos file
        colFermi:   column of Fermi energy starting from the last col   
        spin:       does the Dos take into account spin
                        - 1: No, divide by 1
                        - 2: Yes, divide by 2
        nheader:    number of header lines in Dos file
        nfooter:    number of footer lines in Dos file
        unitFlag:   units used in Dos file  
                        - 0: meV
                        - 1: eV  

    --------------------------------------------------------------------
    Output:
        energies:   Nintx1 vector containing energy grid points, [meV]
        dos:        Nintx1 vector containing dos at epsilon grid points,
                    [1/meV]
        ef:         Fermi energy
        idx_ef:     fermi energy index in energies
    
    --------------------------------------------------------------------
    Comments:
        - Currently the energies must be in column 1 and the dos in 
          column 2 !!
    --------------------------------------------------------------------
    """


    ### Read in dos file ###
    dos_data = readdlm(dos_file) 

    ### Default values ###
    (nheader != -1) || (nheader = findfirst(isa.(dos_data[:, 1], Number)) - 1)
    (nfooter != -1) || (nfooter = size(dos_data, 1) - findlast(isa.(dos_data[:, 1], Number)))

    ### Remove header & footer
    header = join(dos_data[1:nheader, :], " ")
    dos = Float64.(dos_data[nheader+1:end-nfooter, 1:2])
    (~isempty(unit)) || (unit = getUnit(header, "Dos"))

    ### extract energies and dos
    energies = dos[:, 1]
    dos = dos[:, 2]
    dos[dos.<0.0] .= 0.0 # set negative dos to 0

    ### spin
    dos = dos / spin

    ### Fermi energy
    if ef == -1
        try
            ef = Float64.((dos_data[1, end-colFermi]))

            ### Convert units ###
            if unit == "meV"
                ef = ef
                energies = energies
                dos = dos
            elseif unit == "eV"
                ef = ef .* 1000
                energies = energies .* 1000
                dos = dos .* 0.001
            elseif unit == "THz"
                ef = ef * THz2meV
                energies = energies .* THz2meV
                dos = dos ./ THz2meV
            elseif unit == "Ry"
                ef = ef .* Ry2meV
                energies = energies .* Ry2meV
                dos = dos ./ Ry2meV
            else
                error("Invalid Unit! Please check the header of the Dos-file and try again!")
            end

        catch ex

            # log file
            print(log_file, "\nERROR while reading the fermi energy from the Dos-file:")
            showerror(log_file, ex)
            print(log_file, "\n\nCheck the header of the Dos-file or enter the fermi energy by hand via ef")
            print(log_file, "\n\nFor further information please refer to the CRASH file\n")
            close(log_file)

            # crash file
            crashFile = open(inp.outdir * "CRASH", "a")
            print(crashFile, "ERROR while reading the fermi energy from the Dos-file:\n")
            print(crashFile, current_exceptions())
            print(crashFile, "\n\nCheck the header of the Dos-file or enter the fermi energy by hand via ef\n\n")
            close(crashFile)

            # Stop
            rethrow(ex)
        end
    else
        ### Convert units ###
        if unit == "meV"
            energies = energies
            dos = dos
        elseif unit == "eV"
            energies = energies .* 1000
            dos = dos .* 0.001
        elseif unit == "THz"
            energies = energies .* THz2meV
            dos = dos ./ THz2meV
        elseif unit == "Ry"
            energies = energies .* Ry2meV
            dos = dos ./ Ry2meV
        else
            error("Invalid Unit! Either set the unit manually via dos_unit or check the header of the Dos-file and try again!")
        end
    end

    ### Shift energies by ef for cDos ###
    energies = energies .- ef

    return energies, dos, ef, unit
end


"""
    readIn_Weep(Weep_file[, unit, nheader, nfooter])

Read in Weep file containing the sreened coulomb interaction.
Weep data must be in column 3
"""
function readIn_Weep(Weep_file, Weep_col=3, unit="", nheader=-1, nfooter=-1)

    ### Read in Weep file ###
    Weep_data = readdlm(Weep_file);

    ### Default values ###
    (nheader != -1) || (nheader = findfirst(isa.(Weep_data[:,1], Number))-1)
    (nfooter != -1) || (nfooter = size(Weep_data,1) - findlast(isa.(Weep_data[:,1], Number)))

    ### Remove header & footer
    header = join(Weep_data[1:nheader,:], " ")
    Weep = Float64.(Weep_data[nheader+1:end-nfooter, Weep_col])
    
    ### reshape to matrix
    Weep = transpose(reshape(Weep, (Int(sqrt(size(Weep, 1))), Int(sqrt(size(Weep, 1))))))

    ### remove outliers from Weep ###
    outlier = findall(Weep .< 0)  	#.|| Weep .> 1e3
    for k in eachindex(outlier)
        Weep[outlier[k]] = 0
    end

    ### Convert ###
    (~isempty(unit)) || (unit = getUnit(header, "Weep"))
    if "meV" == unit    # meV
        Weep = Weep
    elseif  "eV" == unit     # eV
        Weep = Weep.*1000
    elseif  "Ry" == unit      # Ry
        Weep = Weep.*Ry2meV
    elseif  "Ha" == unit     # Hartree
        Weep = Weep.*Ry2meV*2
    else
        error("Invalid Unit! Please checkt the header of the Weep-file and try again!")
    end

    return Weep, unit

end

"""
    readIn_Wen(Wen_file[, ef, colFermi, unit, nheader, nfooter])

Read in Weep file containing the sreened coulomb interaction.
Weep data must be in column 3
"""
function readIn_Wen(Wen_file, ef=-1, colFermi=0, unit="", nheader=-1, nfooter=-1)
    ### Read in Weep file ###
    Wen_data = readdlm(Wen_file);

    ### Default values ###
    (nheader != -1) || (nheader = findfirst(isa.(Wen_data[:,1], Number))-1)
    (nfooter != -1) || (nfooter = size(Wen_data, 1) - findlast(isa.(Wen_data[:, 1], Number)))

    ### Remove header & footer
    header = join(Wen_data[1:nheader, :], " ")
    W_en = Float64.(Wen_data[nheader+1:end-nfooter, 1])
    (~isempty(unit)) || (unit = getUnit(header, "Wen"))

    ### Fermi energy
    if ef == -1
        try
            ef = Float64.((Wen_data[1, end-colFermi]))  

            ### Convert ###
            if "meV" == unit    # meV
                W_en = W_en
            elseif "eV" == unit     # eV
                ef = ef .* 1000
                W_en = W_en .* 1000
            elseif "Ry" == unit      # Ry
                ef = ef .* Ry2meV
                W_en = W_en .* Ry2meV
            elseif "Ha" == unit     # Hartree
                ef = ef .* Ry2meV * 2
                W_en = W_en .* Ry2meV * 2
            else
                error("Invalid Unit! Please checkt the header of the Weep-file and try again!")
            end

        catch ex

            # log file
            print(log_file, "\nERROR while reading the fermi energy from the Wen-file:")
            showerror(log_file, ex)
            print(log_file, "\n\nCheck the header of the Wen-file or enter the fermi energy by hand via ef")
            print(log_file, "\n\nFor further information please refer to the CRASH file\n")
            close(log_file)

            # crash file
            crashFile = open(inp.outdir * "CRASH", "a")
            print(crashFile, "ERROR while reading the fermi energy from the Wen-file:\n")
            print(crashFile, current_exceptions())
            print(crashFile, "\n\nCheck the header of the Wen-file or enter the fermi energy by hand via ef\n\n")
            close(crashFile)

            # Stop
            rethrow(ex)
        end
    else
        ### Convert only Wen, ef already in correct unit
        if "meV" == unit    # meV
            W_en = W_en
        elseif "eV" == unit     # eV
            W_en = W_en .* 1000
        elseif "Ry" == unit      # Ry
            W_en = W_en .* Ry2meV
        elseif "Ha" == unit     # Hartree
            W_en = W_en .* Ry2meV * 2
        else
            error("Invalid Unit! Please checkt the header of the Weep-file and try again!")
        end
    end


    ### Shift energies by ef ###
    W_en = W_en .- ef


    return W_en

end


# Convert units
function getUnit(header, file=nothing)
    """
    Return true if the specified unit exists in the file,
    e.g. check if "eV" occurs in Dos file header

    -------------------------------------------------------------------
    Input:
        unit:       Energy unit given as string      
        header:     header of a given file  

    --------------------------------------------------------------------
    Output:
        unitFlag:   true if correct unit, false otherwise
    
    --------------------------------------------------------------------
    Comments:
        - When using this function it is important to check first composed 
          units like e.g. meV before eV
          If a string contains meV and one checkes for eV the function will
          return true
    --------------------------------------------------------------------
    """

    units = ["meV", "eV", "THz", "Ry", "Ha"]
    for unit in units
        if occursin(unit, header)
            return unit
        end
    end
    
    println("Auto-extraction of unit from "*file*"-file failed! Please type the correct unit case sensitive into the console:")
    unit = readline();

    return unit

    
end


"""
    discardZeros(Dos, energies)

Discard zeros in dos.
"""
function discardZeros(Dos::Vector{Float64}, energies::Vector{Float64})
    idxLower = findfirst(Dos .!= 0)
    idxUpper = findlast(Dos .!= 0)

    return     Dos[idxLower:idxUpper], energies[idxLower:idxUpper]
end



"""
    restrictInput(n, enWndw, ef, Dos, energies[, Weep])

Restrict the input (Weep, Dos and energies) to a energy window around the fermi energy and a subset of grid points

# Examples
```julia-repl
julia> restrictInput(30, [2,3], [10], [pi.*sqrt.(range(5,15,200))], [range(5,15,200)])
```
"""
function restrictInput(n::Int, enWndw::Any, Dos::Vector{Float64}, energies::Vector{Float64}, Weep::Matrix{Float64})

    (~isnothing(enWndw)) || (enWndw = [ -energies[1], energies[end]])
  

    idxLower = findfirst(- enWndw[1] .<= energies)
    idxUpper = findlast(enWndw[2] .>= energies)

    # total grid points
    N = idxUpper - idxLower

    if n > N
        Weep     = Weep[idxLower:idxUpper, idxLower:idxUpper]
        Dos      = Dos[idxLower:idxUpper] 
        energies = energies[idxLower:idxUpper]

    else
        Ninterval = N / n
        rngRed = Int.(round.(collect(idxLower:Ninterval:idxUpper)))

        Weep = Weep[rngRed, rngRed]
        Dos = Dos[rngRed]
        energies = energies[rngRed]

    end

    return Weep, Dos, energies

end


function restrictInput(n::Int, enWndw::Any, Dos::Vector{Float64}, energies::Vector{Float64})
 
    (~isnothing(enWndw)) || (enWndw = [-energies[1], energies[end]])
  

    idxLower = findfirst(-enWndw[1] .<= energies)
    idxUpper = findlast(enWndw[2] .>= energies)

    # total grid points
    N = idxUpper - idxLower

    if n > N
        Dos      = Dos[idxLower:idxUpper] 
        energies = energies[idxLower:idxUpper]

    else
        Ninterval = N / n
        rngRed = Int.(round.(collect(idxLower:Ninterval:idxUpper)))

        Dos = Dos[rngRed]
        energies = energies[rngRed]

    end

    return Dos, energies

end


