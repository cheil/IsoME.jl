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
    InputParser(inp, log_file)

Read, convert, preprocess inputs.
"""
function InputParser(inp::arguments, log_file)

    ### Init table size ###
    console = Dict()
    if inp.include_Weep == 1 && inp.cDOS_flag == 0
        # header table
        console["header"] = ["it", "phic", "phiph", "znormi", "shifti", "ef-mu", "deltai", "err_delta"]
        # width table
        console["width"] = [8, 10, 10, 10, 10, 10, 10, 11]
        # precision data
        console["precision"] = [0, 2, 2, 2, 2, 2, 2, 5]

    elseif inp.include_Weep == 1 && inp.cDOS_flag == 1
        # header table
        console["header"] = ["it", "phic", "phiph", "znormi", "deltai", "err_delta"]
        #width table
        console["width"] = [8, 10, 10, 10, 10, 11]
        # precision data
        console["precision"] = [0, 2, 2, 2, 2, 5]

    elseif inp.include_Weep == 0 && inp.cDOS_flag == 0
        # header table
        console["header"] = ["it", "znormi", "shifti", "ef-mu", "deltai", "err_delta"]
        #width table
        console["width"] = [8, 10, 10, 11, 10, 11]
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
        error("Unkwon mode! Check if the cDOS_flag and include_Weep flag are set correctly!")
    end
    console = formatTableHeader(console)

    console = printStartMessage(console, log_file)
  

    ########## READ-IN ##########
    ####### a2f file #######
    a2f_omega, a2f, inp.ind_smear, inp.a2f_unit = readIn_a2f(inp.a2f_file, inp.ind_smear, inp.a2f_unit, inp.nheader_a2f, inp.nfooter_a2f, inp.nsmear)


    ########## Dos and Weep ##########
    if isfile(inp.dos_file)     # ~isempty was wrong
        # read dos
        dos_en, dos, ef, inp.dos_unit = readIn_Dos(inp.dos_file, inp.ef, inp.spinDos, inp.dos_unit, inp.nheader_dos, inp.nfooter_dos, outdir = inp.outdir, logFile = log_file)
        inp.ef = ef

        # remove zeros at begining/end of dos
        dos, dos_en = discardZeros(dos, dos_en)

        ### Interpolation ###
        # Interpolation Object DoS
        epsilonItp = range(dos_en[1], dos_en[end], length(dos_en))      
        itpDos = scale(interpolate(dos, BSpline(Cubic())), epsilonItp) 

        # encut in meV
        if (inp.encut > dos_en[end]) || (-inp.encut < dos_en[1])    
            text = "Energy cutoff exceeds range of DOS!"
            printWarning(text, log_file)
        end

        if inp.include_Weep == 1 || ((isfile(inp.Weep_file) && (isempty(inp.Wen_file) || isfile(inp.Wen_file))) && inp.mu == -1)  
            # read Weep + energy grid points
            Weep, Wen, inp.efW, inp.Weep_unit = readIn_Weep(inp.Weep_file, inp.Wen_file, inp.Weep_col, inp.Wen_col, inp.efW, inp.Weep_unit, inp.nheader_Weep, inp.nfooter_Weep, inp.nheader_Wen, inp.nfooter_Wen, outdir = inp.outdir, logFile = log_file)

            ### Interpolation ###
            # Interpolation Object Weep
            epsilonItp = range(Wen[1], Wen[end], length(Wen))       
            itpWeep = scale(interpolate(Weep, BSpline(Linear())), (epsilonItp, epsilonItp)) 

            # interpolate 
            dos_en, dos, Weep = interpolateInputs(itpDos, dos_en, inp.itpStepSize, inp.itpBounds, inp.encut, itpWeep = itpWeep, Wen = Wen)

            # mu, idxWef = idx_ef
            (inp.mu != -1) || (idxWef = findmin(abs.(dos_en))[2]; inp.mu = dos[idxWef].*Weep[idxWef,idxWef])

        else
            # interpolate 
            dos_en, dos, Weep = interpolateInputs(itpDos, dos_en, inp.itpStepSize, inp.itpBounds, inp.encut)
        end

        # idx 
        idxShiftcut = [findfirst(dos_en .> -inp.shiftcut), findlast(dos_en .< inp.shiftcut)]

    else
        # default values
        dos = []
        dos_en = []
        ef = nothing
        idxShiftcut = -1
        Weep = nothing
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

    ### calc mu*'s
    if inp.mu == -1 && inp.muc_ME == -1 && inp.muc_AD == -1
        # defaut value
        inp.muc_AD = 0.12
        calcMucME(inp, a2f, a2f_omega, log_file)

    elseif inp.muc_AD == -1 && inp.muc_ME == -1
        if inp.typEl != -1
            calcMucs(inp, inp.typEl, a2f, a2f_omega, log_file)
        
        elseif ~(isnothing(ef) || ef == -1 || ef == 0)
            inp.typEl = ef
            calcMucs(inp, inp.typEl, a2f, a2f_omega, log_file)

        elseif ~(inp.efW == -1 || inp.efW == 0)
            inp.typEl = inp.efW
            calcMucs(inp, inp.typEl, a2f, a2f_omega, log_file)

        else
            text = "Unable to calculate μ* from μ without a typical electron energy!"
            text *= "\nConsider setting typEl, ef or efW manually."
            text *= "\nUsing μ*_AD = 0.12 instead"
            printWarning(text, log_file)

            inp.muc_AD = 0.12
            calcMucME(inp, a2f, a2f_omega, log_file)
        end

    elseif  inp.include_Weep == 0 && inp.muc_AD != -1 && inp.muc_ME == -1   
        calcMucME(inp, a2f, a2f_omega, log_file)

    elseif  inp.muc_ME != -1 && inp.muc_AD == -1 # Wenn nur μ*_ME gegeben wird kein μ*_AD berechnet falls mit W und μ* sollen nur berechnet werden falls -1
        calcMucAD(inp, a2f, a2f_omega)
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
    matval = (a2f_omega, a2f, dos_en, dos, Weep, dosef, idx_ef, ndos, BCS_gap, idxShiftcut)

    return inp, console, matval, ML_Tc
end


"""
    createDirectory(inp, strIsoME)

Create directory.
"""
function createDirectory(inp::arguments, strIsoME::String)

    if inp.testMode # hidden outputs in test mode
        log_file = IOBuffer()
        errorLogger = SimpleLogger(log_file, Logging.Error)
    else
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
            mkpath(inp.outdir)   # mkpath
        catch ex
            error("Couldn't write into " * inp.outdir * "! Outdir may not writable or an invalid path.\n\n")
        end

        log_file = open(inp.outdir * "log.txt", "w")
        print(log_file, strIsoME)

        # logging to console and log-file (@warn,...)
        errorLogger = SimpleLogger(log_file, Logging.Error)
        file_logger = SimpleLogger(log_file)
        tee_logger = TeeLogger(ConsoleLogger(), file_logger)
        global_logger(tee_logger)
    end

    return inp, log_file, errorLogger

end


"""
    checkInput!(inp)

Check which input files (a2f, dos, weep) exist.
"""
function checkInput(inp::arguments)

    # check input files / cDOS & Weep
    if ~isfile(inp.a2f_file)
        error("Invalid path to a2f-file!")
        
    elseif ~isfile(inp.dos_file) && (inp.cDOS_flag == 0 || inp.include_Weep == 1)

        text = "Invalid path to Dos-file!\n\n"
        error(text)
    
    elseif inp.include_Weep == 1 && ((~isfile(inp.Weep_file)) || (~isempty(inp.Wen_file) && ~isfile(inp.Wen_file)))
        
        text = "Invalid path to Weep or Wen-file!\n\n"
        error(text)

    end

    return inp

end

"""
    readIn_a2f(a2f_file, indSmear=-1, unit="", nheader=-1, nfooter=-1, nsmear=-1)   

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
    elseif  "Ha" == unit     # Hartree
        omega_raw = omega_raw .* Ry2meV*2
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


    return omega_fine, a2f_fine, indSmear, unit

end


# Read in DOS
"""
    readIn_Dos(dos_file, ef =-1, spin=2, unit="", nheader=-1, nfooter=-1; outdir ="./", logFile = nothing)

Read the dos file. All quantities are converted to meV.

The energies must be in column 1 and the dos in column 2
"""
function readIn_Dos(dos_file, ef =-1, spin=2, unit="", nheader=-1, nfooter=-1; outdir ="./", logFile = nothing)

    ### Read in dos file ###
    dos_data = readdlm(dos_file) 

    ### Default values ###
    (nheader != -1) || (nheader = findfirst(isa.(dos_data[:, 1], Number)) - 1)
    (nfooter != -1) || (nfooter = size(dos_data, 1) - findlast(isa.(dos_data[:, 1], Number)))

    ### Remove header & footer
    header = dos_data[1:nheader, :]
    dos = Float64.(dos_data[nheader+1:end-nfooter, 1:2])
    (~isempty(unit)) || (unit = getUnit(join(header, " "), "Dos"))

    ### extract energies and dos
    energies = dos[:, 1]
    dos = dos[:, 2]
    dos[dos.<0.0] .= 0.0 # set negative dos to 0

    ### spin
    dos = dos / spin

    ### Fermi energy
    if ef == -1
        ef = extractFermiEnergy(header, unit, "Weep", outdir = outdir, logFile = logFile)
    end

    ### Convert 
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
    elseif  "Ha" == unit     # Hartree
        energies = Weep.*Ry2meV*2
        dos = Wen ./ Ry2meV*2
    else
        error("Invalid Unit! Either set the unit manually via dos_unit or check the header of the Dos-file and try again!")
    end

    ### Shift energies by ef for cDos ###
    energies = energies .- ef

    return energies, dos, ef, unit
end


"""
    readIn_Weep(Weep_file, Wen_file="", Weep_col=3, Wen_col=1, ef=-1, unit = "", nheader=-1, nfooter=-1,  nheaderWen=-1, nfooterWen=-1; outdir = "./", logFile = nothing)

Read in Weep file containing the sreened coulomb interaction.
Weep data must be in column 3
"""
function readIn_Weep(Weep_file, Wen_file="", Weep_col=3, Wen_col=1, ef=-1, unit = "", nheader=-1, nfooter=-1,  nheaderWen=-1, nfooterWen=-1; outdir = "./", logFile = nothing)

    ### Read in Weep file ###
    Weep_data = readdlm(Weep_file);

    # Default values
    (nheader != -1) || (nheader = findfirst(isa.(Weep_data[:,1], Number))-1)
    (nfooter != -1) || (nfooter = size(Weep_data,1) - findlast(isa.(Weep_data[:,1], Number)))

    # Remove header & footer
    header = Weep_data[1:nheader,:]
    Weep = Float64.(Weep_data[nheader+1:end-nfooter, Weep_col])
    
    # reshape to matrix
    numWens = Int(sqrt(size(Weep, 1)))
    Weep = transpose(reshape(Weep, numWens, numWens))

    # remove outliers from Weep 
    Weep[Weep .< 0] .= 0

    # unit
    (~isempty(unit)) || (unit = getUnit(join(header, " "), "Weep"))

    ### Fermi energy
    if ef == -1
        ef = extractFermiEnergy(header, unit, "Weep", outdir = outdir, logFile = logFile)
    end

    ### read in W energies ###
    if isempty(Wen_file)
        #Wen = Float64.(Weep_data[nheader+1:numWens+nheader, Wen_col])
        Wen = unique(Float64.(Weep_data[nheader+1:end-nfooter, Wen_col]))
    else
        Wen = readIn_Wen(Wen_file, Wen_col, nheaderWen, nfooterWen)
    end

    ### Convert ###
    if "meV" == unit    # meV
        Weep = Weep
        Wen = Wen
    elseif  "eV" == unit     # eV
        Weep = Weep.*1000
        Wen = Wen.*1000
    elseif  "Ry" == unit      # Ry
        Weep = Weep.*Ry2meV
        Wen = Wen.*Ry2meV
    elseif  "Ha" == unit     # Hartree
        Weep = Weep.*Ry2meV*2
        Wen = Wen.*Ry2meV*2
    elseif unit == "THz"
        Weep = Weep .* THz2meV
        Wen = Wen .* THz2meV
    else
        error("Invalid Unit! Consider setting the unit manually (Weep_unit) or check the header of the Weep-file and try again!")
    end

    # shift by ef
    Wen = Wen .- ef

    return Weep, Wen, ef, unit

end

"""
    readIn_Wen(Wen_file, Wen_col, nheader=-1, nfooter=-1)

Read in Wen 
Energy grid for Weep 
"""
function readIn_Wen(Wen_file, Wen_col, nheader=-1, nfooter=-1)
    ### Read in Weep file ###
    Wen_data = readdlm(Wen_file);

    ### Default values ###
    (nheader != -1) || (nheader = findfirst(isa.(Wen_data[:,1], Number))-1)
    (nfooter != -1) || (nfooter = size(Wen_data, 1) - findlast(isa.(Wen_data[:, 1], Number)))

    ### Remove header & footer
    Wen = Float64.(Wen_data[nheader+1:end-nfooter, Wen_col])

    return Wen

end


# Fermi energy
"""
    extractFermiEnergy(header, unit, nameFile=nothing)

Extract the fermi energy from the header of the input files
"""
function extractFermiEnergy(header, unit, nameFile=nothing; outdir = "./", logFile = nothing)

    ef = -1
    try
        logNums = isa.(header, Number)
        if sum(logNums) == 1 
            ef = Float64(only(header[logNums]))
        elseif sum(logNums[:, 2:end] .& (header[:, 1:end-1] .== "=" )) == 1 
            ef = Float64(only(header[:, 2:end][logNums[:, 2:end] .& (header[:, 1:end-1] .== "=" )]))
        else
            nameFermi = ["efermi", "ef", "fermi", "e_fermi"]
            lcHeader = map(x -> isa(x, AbstractString) ? lowercase(x) : x, header)

            for name in nameFermi
                if any(lcHeader .== name)
                    idx = findall(lcHeader .== name)
                    row = idx[1][1]
                    col = idx[1][2]
                    ef = Float64.(only(header[row, findfirst(isa.(header[row, col:end], Number))+col-1]))
                    break
                end
            end 
        end

        ### Convert ###
        if "meV" == unit # meV
            ef = ef
        elseif "eV" == unit     # eV
            ef = ef .* 1000
        elseif "Ry" == unit      # Ry
            ef = ef .* Ry2meV
        elseif "Ha" == unit     # Hartree
            ef = ef .* Ry2meV * 2
        else
            error("Invalid Unit! Consider setting the unit manually ("*nameFile*"_unit) or check the header of the "*nameFile*"-file and try again!")
        end

    catch ex
        text = "Error while reading the fermi energy from the " * nameFile * "-file."
        text *= "\nConsider setting the fermi-energy manually (ef or efW) or check the header of the " * nameFile * "-file\n\n"
        @error text
    end

    return ef
end


# Convert units
"""
    getUnit(header, nameFile=nothing)

Read the units from the header of the input files
"""
function getUnit(header, nameFile=nothing)

    units = ["meV", "eV", "THz", "Ry", "Ha"]
    for unit in units
        if occursin(unit, header)
            return unit
        end
    end
    
    println("Auto-extraction of unit from "*nameFile*"-file failed! Please type the correct unit case sensitive into the console (it might be that you need to type it twice due to a bug in julia):")
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
    calcMucME(inp, console, a2f_omega, log_file)

Calculate μ*_ME from μ*_AD using formula as in 
Pellegrini, Ab initio methods for superconductivity 
DOI: 10.1038/s42254-024-00738-9
"""
function calcMucME(inp, a2f, a2f_omega, log_file)
    inp.muc_ME = inp.muc_AD / (1 + inp.muc_AD * log(maximum(a2f_omega[a2f.>0.01]) / inp.omega_c))

    if inp.muc_ME < 0 || inp.muc_ME > 0.8 || inp.muc_ME > 3 * inp.muc_AD
        inp.muc_ME = minimum([3 * inp.muc_AD, 0.8])

        text = "Couldn't calculate a reasonable μ*_ME from μ*_AD."
        text *= "\nUsing μ*_ME = minimum(3*μ*_AD, 0.8) instead."
        text *= "\nConsider setting it manually!"
        printWarning(text, log_file)
    end
end


"""
    calcMucAD(inp, console, a2f_omega)

Calculate μ*_AD from μ*_ME using formula (30) in
Pellegrini, Ab initio methods for superconductivity 
DOI: 10.1038/s42254-024-00738-9
"""
function calcMucAD(inp, a2f, a2f_omega)
    inp.muc_AD = inp.muc_ME / (1 - inp.muc_ME * log(maximum(a2f_omega[a2f.>0.01]) / inp.omega_c))

    if inp.muc_AD < 0 || inp.muc_AD > 0.2 
        inp.muc_AD = 0.12   # default
    end
end


"""
    calcMucME(inp, console, a2f_omega, log_file)

Calculate μ*_ME and μ*_AD from μ using formulas as in 
Pellegrini, Ab initio methods for superconductivity 
DOI: 10.1038/s42254-024-00738-9
"""
function calcMucs(inp, ef, a2f, a2f_omega, log_file)
    # μ*_ME < 4*μ
    if inp.omega_c > ef * exp(3 / (4 * inp.mu))
        inp.omega_c = ef * exp(3 / (4 * inp.mu))

        text = "Matsubara cutoff would lead to μ*_ME > 4*μ."
        text *= "\nomega_c has been set to a smaller value."
        text *= "\nConsider checking the typical electronic energy typEl!"
        printWarning(text, log_file)
    end

    inp.muc_AD = inp.mu / (1 + inp.mu * log(ef / maximum(a2f_omega[a2f.>0.01])))
    if inp.include_Weep == 0
        inp.muc_ME = inp.mu / (1 + inp.mu * log(ef / inp.omega_c))
    end
end
