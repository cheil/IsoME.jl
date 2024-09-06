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
function InputParser(inp::arguments)

    inp = checkInput!(inp)

    ### open log_file ###
    if inp.flag_log == 1
        if isfile(inp.outdir * "log.txt")
            rm(inp.outdir * "log.txt")
        end

        inp.log_file = open(inp.outdir * "log.txt", "w")
    end

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
        console["width"] = [8, 10, 10, 10, 10, 12]
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

    console = printStartMessage(inp, console)
  

    ########## READ-IN ##########
    ####### a2f file #######
    a2f_omega_fine, a2f_fine = readIn_a2f(inp.a2f_file, inp.ind_smear, inp.a2f_unit, inp.nheader_a2f, inp.nfooter_a2f, inp.nsmear)

    # determine superconducting properties from Allen-Dynes McMillan equation based on interpolated a2F
    AD_data = calc_AD_Tc(a2f_omega_fine, a2f_fine, inp.muc_AD)
    ML_Tc = AD_data[1] / kb;    # ML-Tc in K
    AD_Tc = AD_data[2] / kb;    # AD-Tc in K
    BCS_gap = AD_data[3];       # BCS gap value in meV
    lambda = AD_data[4];        # total lambda
    omega_log = AD_data[5];     # omega_log in meV


    ########## read-in and process Dos and Weep ##########
    if ~isempty(inp.dos_file) 
         # read dos
         dos_en, dos, ef, inp.dos_unit = readIn_Dos(inp.dos_file, inp.cDOS_flag, inp.colFermi_dos, inp.spinDos, inp.dos_unit, inp.nheader_dos, inp.nfooter_dos)
            
        if inp.include_Weep == 1
            # read Weep + energy grid points
            Weep, inp.Weep_unit = readIn_Weep(inp.Weep_file, inp.Weep_col, inp.Weep_unit, inp.nheader_Weep, inp.nfooter_Weep)
            W_en = readIn_Wen(inp.Wen_file, inp.Wen_unit, inp.nheader_Wen, inp.nfooter_Wen)

            # overlap of energies
            en_interval = [W_en[findfirst(W_en .> dos_en[1])]; W_en[findlast(W_en .< dos_en[end])]]
            # number of points overlapping
            idxOverlap = [findfirst(W_en .> dos_en[1]), findlast(W_en .< dos_en[end])]
            Nitp = idxOverlap[2] - idxOverlap[1] + 1
            
            # restrict Weep
            Weep = Weep[idxOverlap[1]:idxOverlap[2], idxOverlap[1]:idxOverlap[2]]
            
            # Interpolation Dos
            dos_en, dos = interpolateDos(dos_en, dos, en_interval, Nitp)
        else
            # no Weep
            Weep = nothing
            W_en = nothing

            # Interpolation Dos
            #dos_en, dos = interpolateDos(dos_en, dos, [dos_en[1], dos_en[end]], 10*length(dos_en))
        end
    else
        dos = []
        dos_en = []
        ef = nothing
    end

    # user specified ef
    if inp.ef != -1
        ef = inp.ef 
    end

    ### calc mu*'s
    if inp.include_Weep == 0
        if inp.mu != -1
            if isnothing(ef)
                print(@yellow "WARNING: ")
                print("Unable to calculate μ* from μ without the fermi-energy! Either provide a dos-file or set the fermi-energy directly in the input structure.\n\n")
            
                if inp.flag_log == 1
                    print(inp.log_file, "WARNING: Unable to calculate μ* from μ without the fermi-energy! Either provide a dos-file or set the fermi-energy directly in the input structure.\n\n")
                end
            else
                # ensure μ*_ME < 4*μ --> reasonable ??
                if inp.omega_c > ef*exp(3/(4*inp.mu))
                    inp.omega_c = ef*exp(3/(4*inp.mu))

                    print(@yellow "WARNING: ")
                    print("Matsubara cutoff would lead to μ*_ME > 4*μ. omega_c has been set to a smaller value.\n\n")
            
                    if inp.flag_log == 1
                        print(inp.log_file, "WARNING: Matsubara cutoff would lead to μ*_ME > 4*μ. omega_c has been set to a smaller value.\n\n")
                    end
                end

                inp.muc_AD = inp.mu / (1 + inp.mu*log(ef/maximum(a2f_omega_fine[a2f_fine .> 0.01])))
                inp.muc_ME = inp.mu / (1 + inp.mu*log(ef/inp.omega_c))
            end
        elseif inp.muc_ME == -1
            inp.muc_ME = inp.muc_AD / (1 + inp.muc_AD*log(maximum(a2f_omega_fine[a2f_fine .> 0.01])/inp.omega_c))

            if inp.muc_ME < 0 || inp.muc_ME > 3*inp.muc_AD || inp.muc_ME > 0.8
                inp.muc_ME = minimum([3*inp.muc_AD, 0.8])

                print(@yellow "WARNING: ")
                print("Couldn't calculate a reasonable μ*_ME from μ*_AD. Using μ*_ME = minimum(3*μ*_AD, 0.8) instead. Consider setting it manually! \n\n")
            
                if inp.flag_log == 1
                    print(inp.log_file, "WARNING: Couldn't find a reasonable μ*_ME from μ*_AD. Using μ*_ME = minimum(3*μ*_AD, 0.8) instead. Consider setting it manually! \n\n")
                end
            end
        end
    end

    # print Allen-Dynes
    printADtable(inp, console, ML_Tc, AD_Tc, BCS_gap, lambda, omega_log)


    # restrict Dos/Weep to a subset of grid points
    if inp.Nrestrict != -1
        if inp.include_Weep == 1
            # call restrict function
            Weep, dos, dos_en = restrictInput(inp.Nrestrict, inp.wndRestrict, ef, dos, dos_en, Weep)
        else
            # call restrict function
            dos, dos_en = restrictInput(inp.Nrestrict, inp.wndRestrict, ef, dos, dos_en)
        end
    end


    ### remove zeros at begining/end of dos ###
    if inp.include_Weep == 1 
        dos, dos_en, Weep = neglectZeros(dos, dos_en, Weep = Weep)
    elseif inp.cDOS_flag == 0
        dos, dos_en = neglectZeros(dos, dos_en)
    end


    if inp.include_Weep == 0 && inp.cDOS_flag == 1
        # default values in case no dos-file is given
        ndos = 0
        idx_ef = 0
        dosef = 0
    else
        # length energy vector
        ndos = size(dos_en, 1)

        # index of fermi energy
        idx_ef = findmin(abs.(dos_en .- ef))
        idx_ef = idx_ef[2]

        # dos at ef
        dosef = dos[idx_ef]
    end

    # material specific values
    matval = (a2f_omega_fine, a2f_fine, dos_en, dos, Weep, ef, dosef, idx_ef, ndos, BCS_gap)

    return inp, console, matval, ML_Tc
end


"""
    checkInput!()

Check which input files (a2f, dos, weep) exist. Overwrite flags if neccessary.
"""
function checkInput!(inp::arguments)
    # input files / cDOS & Weep
    if ~isfile(inp.a2f_file)
        error("Invalid path to a2f-file!")
        
    elseif ~isfile(inp.dos_file) && (inp.cDOS_flag == 0 || inp.include_Weep == 1)
        inp.cDOS_flag = 1 
        inp.include_Weep = 0
        print(@yellow "WARNING: ")
        print("No valid Dos-file specified! Calculating Tc within constant Dos approximation using Anderson-pseudopotential instead\n\n")
    
        if inp.flag_log == 1
            print(inp.log_file, "WARNING: No valid Dos-file specified! Calculating Tc within constant Dos approximation using Anderson-pseudopotential instead\n\n")
        end
    
    elseif (~(isfile(inp.Weep_file) || ~isfile(inp.Wen_file))) && inp.include_Weep == 1
        inp.include_Weep = 0
        print(@yellow "WARNING: ")
        print("No valid Weep/Energies-file specified! Calculating Tc using Anderson-pseudopotential instead\n\n")
    
        if inp.flag_log == 1
            print(inp.log_file, "WARNING: No valid Weep/DosWeep-file specified! Calculating Tc using Anderson-pseudopotential instead\n\n")
        end
    end


    # Tc search mode
    if inp.temps == [-1]
        inp.TcSearchMode_flag = 1
        print(@yellow "WARNING: ")
        print("No temperatures specified - Activating automatic Tc search mode instead!\n\n")
    
        if inp.flag_log == 1
            print(inp.log_file, "WARNING: No temperatures specified - Activating automatic Tc search mode instead!\n\n")
        end
    end


    # create output directory
    if ~(inp.outdir[end] == '/' || inp.outdir[end] == '\\')
        inp.outdir = inp.outdir*"/"
    end

    if ~isdir(inp.outdir)
        try
            mkdir(inp.outdir)
        catch
            error("Couldn't write into "*inp.outdir*"! Path in outdir may not be writable or invalid.")
        end
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
    (nsmear != -1) || (nsmear = length(a2f_data[nheader, isa.(a2f_data[nheader,:], Number)])-1)
    (indSmear != -1) || (indSmear = Int64(ceil(nsmear/2)))

    ### Remove header & footer
    header = join(a2f_data[1:nheader,:], " ")
    a2f_data = Float64.(a2f_data[nheader+1:end-nfooter, 1:nsmear+1])  # previous version: nheader+1:end-nfooter

    ### Convert omega ###
    omega_raw = a2f_data[:, 1]
    (~isempty(unit)) || (unit = getUnit(header, "a2F"))
    if "meV" == unit
        omega_raw = omega_raw
    elseif "THz" == unit
            omega_raw = omega_raw * THz2meV
    elseif "Ry" == unit
        omega_raw = omega_raw * Ry2meV
    else
        error("Invalid Unit! Please checkt the header of the a2F-file and try again!")
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
function readIn_Dos(dos_file, cDOS_flag, colFermi=1, spin=2, unit="", nheader=-1, nfooter=-1)
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
    (nheader != -1) || (nheader = findfirst(isa.(dos_data[:,1], Number))-1)
    (nfooter != -1) || (nfooter = size(dos_data,1) - findlast(isa.(dos_data[:,1], Number)))

    ### Fermi energy
    ef = Float64.((dos_data[1, end-colFermi]))  

    ### Remove header & footer
    header = join(dos_data[1:nheader,:], " ")
    dos_data = Float64.(dos_data[nheader+1:end-nfooter, 1:2])

    ### extract energies and dos
    energies = dos_data[:,1]
    dos = dos_data[:,2]
    dos[dos.<0.0] .= 0.0 # set negative dos to 0

    ### spin
    dos = dos/spin

    ### Convert units ###
    (~isempty(unit)) || (unit = getUnit(header, "Dos"))
    if unit == "meV"
        ef = ef
        energies = energies
        dos = dos
    elseif unit == "eV"
        ef = ef .*1000 
        energies = energies .*1000
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


    ### Shift energies by ef for cDos ###
    if cDOS_flag == 1
        energies = energies .- ef
    end

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
    outlier = findall(Weep .< 0)    #.|| Weep .> 1e5	
    for k in eachindex(outlier)
        Weep[outlier[k]] = 0
    end

    ### Convert ###
    (~isempty(unit)) || (unit = getUnit(header, "Weep"))
    if "meV" == unit    # meV
        Weep = Weep
    elseif  "eV" == unit     # eV
        Weep = Weep.*1000
    else
        error("Invalid Unit! Please checkt the header of the Weep-file and try again!")
    end

    ### make Weep symmetric ###
    # Weep should be symmetric anyway
    #Weep = Symmetric(Weep)

    return Weep, unit

end

"""
    readIn_Weep(Weep_file[, unit, nheader, nfooter])

Read in Weep file containing the sreened coulomb interaction.
Weep data must be in column 3
"""
function readIn_Wen(Wen_file, unit="", nheader=-1, nfooter=-1)
    ### Read in Weep file ###
    Wen_data = readdlm(Wen_file);

    ### Default values ###
    (nheader != -1) || (nheader = findfirst(isa.(Wen_data[:,1], Number))-1)
    (nfooter != -1) || (nfooter = size(Wen_data,1) - findlast(isa.(Wen_data[:,1], Number)))

    ### Remove header & footer
    header = join(Wen_data[1:nheader,:], " ")
    W_en = Float64.(Wen_data[nheader+1:end-nfooter, 1])

    ### Convert ###
    (~isempty(unit)) || (unit = getUnit(header, "Wen"))
    if "meV" == unit    # meV
        W_en = W_en
    elseif  "eV" == unit     # eV
        W_en = W_en.*1000
    else
        error("Invalid Unit! Please checkt the header of the Weep-file and try again!")
    end

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

    units = ["meV", "eV", "THz", "Ry"]
    for unit in units
        if occursin(unit, header)
            return unit
        end
    end
    
    println("Auto-extraction of unit from "*file*"-file failed! Please type the correct unit TWICE(!) into the console:")
    unit = readline();

    return unit

    
end


"""
    neglectZeros(Dos, energies[, Weep])

Discard zeros in dos.
"""
function neglectZeros(Dos::Vector{Float64}, energies::Vector{Float64}; Weep = nothing)
    idxLower = findfirst(Dos .!= 0)
    idxUpper = findlast(Dos .!= 0)

    if isnothing(Weep)
        return     Dos[idxLower:idxUpper], energies[idxLower:idxUpper]
    else
        return     Dos[idxLower:idxUpper], energies[idxLower:idxUpper], Weep[idxLower:idxUpper, idxLower:idxUpper]
    end
end



"""
    restrictInput(n, enWndw, ef, Dos, energies[, Weep])

Restrict the input (Weep, Dos and energies) to a energy window around the fermi energy and a subset of grid points

# Examples
```julia-repl
julia> restrictInput(30, [2,3], [10], [pi.*sqrt.(range(5,15,200))], [range(5,15,200)])
```
"""
function restrictInput(n::Int, enWndw::Any, ef::Float64, Dos::Vector{Float64}, energies::Vector{Float64}, Weep::Matrix{Float64})

    (~isnothing(enWndw)) || (enWndw = [ef - energies[1], energies[end] - ef])
  

    idxLower = findfirst(ef - enWndw[1] .<= energies)
    idxUpper = findlast(ef + enWndw[2] .>= energies)

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


function restrictInput(n::Int, enWndw::Any, ef::Float64, Dos::Vector{Float64}, energies::Vector{Float64})
 
    (~isnothing(enWndw)) || (enWndw = [ef - energies[1], energies[end] - ef])
  

    idxLower = findfirst(ef - enWndw[1] .<= energies)
    idxUpper = findlast(ef + enWndw[2] .>= energies)

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


