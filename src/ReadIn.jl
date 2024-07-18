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

#=
struct iso
    # Parameters
    muc = muc
    temps = temps
    omega_c = omega_c

    # mode
    cDOS_flag::Int64 = cDOS_flag    # boolean?
    TcSearchMode_flag::Int64 = TcSearchMode_flag
    mu_flag::Int64 = mu_flag

    # a2f input file
    ind_smear::Int64    = ind_smear
    nsmear::Int64       = nsmear
    nheader_a2f::Int64  = nheader_a2f
    nfooter_a2f::Int64  = nfooter_a2f
    a2f_unit::Int64     = a2f_unit

    # dos input file
    nheader_dos::Int64  = nheader_dos 
    nfooter_dos::Int64  = nfooter_dos 
    dos_unit::Int64     = dos_unit 
    spinDos::Int64      = spinDos 
    colFermi_dos::Int64 = colFermi_dos 
    nheader_dosW::Int64 = nheader_dosW 
    nfooter_dosW::Int64 = nfooter_dosW 
    dosW_unit::Int64    = dosW_unit 
    spinDosW::Int64     = spinDosW 
    colFermi_dosW::Int64 = colFermi_dosW 

    # Weep input file, auto extraction if nothing
    nheader_Weep::Int64 = nheader_Weep 
    nfooter_Weep::Int64 = nfooter_Weep 
    Weep_unit::Int64 = Weep_unit 
    include_Weep::Int64 = include_Weep

    # Output
    outdir::string = outdir
    flag_log::Int64 = flag_log
    flag_figure::Int64 = flag_figure 
    flag_outfile::Int64 = flag_outfile 

    # Restrict Weep, REMOVE BEFORE MERGE !!!
    Nrestrict::Int64 = Nrestrict 
    wndRestrict::Vector{Float64} = wndRestrict 


end
=#

"""
    InputParser()

Include the file specified via path.

# Examples
"""
function InputParser(path)
    include(path)
        
    ### Defining some defaults ###s
    # Parameters
    (@isdefined muc) || (muc = 0.1)
    (@isdefined temps) || (temps = collect(5.0:5.0:200.0))
    (@isdefined omega_c) || (omega_c = 20.0 * a2f_omega_raw[end])

    # mode
    (@isdefined cDOS_flag) || (cDOS_flag = 0)
    (@isdefined TcSearchMode_flag) || (TcSearchMode_flag = 0)
    (@isdefined mu_flag) || (mu_flag = 1)

    # a2f input file, auto extraction if nothing
    (@isdefined ind_smear) || (ind_smear = 1)
    (@isdefined nsmear) || (nsmear = nothing)
    (@isdefined nheader_a2f) || (nheader_a2f = nothing)
    (@isdefined nfooter_a2f) || (nfooter_a2f = nothing)
    (@isdefined a2f_unit) || (a2f_unit = nothing)

    # dos input file, auto extraction if nothing
    (@isdefined nheader_dos) || (nheader_dos = nothing)
    (@isdefined nfooter_dos) || (nfooter_dos = nothing)
    (@isdefined dos_unit) || (dos_unit = nothing)
    (@isdefined spinDos) || (spinDos = 2)
    (@isdefined colFermi_dos) || (colFermi_dos = 1)
    (@isdefined nheader_dosW) || (nheader_dosW = nothing)
    (@isdefined nfooter_dosW) || (nfooter_dosW = nothing)
    (@isdefined dosW_unit) || (dosW_unit = nothing)
    (@isdefined spinDosW) || (spinDosW = 1)
    (@isdefined colFermi_dosW) || (colFermi_dosW = 0)

    # Weep input file, auto extraction if nothing
    (@isdefined nheader_Weep) || (nheader_Weep = nothing)
    (@isdefined nfooter_Weep) || (nfooter_Weep = nothing)
    (@isdefined Weep_unit) || (Weep_unit = nothing)
    (@isdefined include_Weep) || (include_Weep = 1)

    # Output
    (@isdefined outdir) || (outdir = "./output/" )
    (@isdefined flag_log) || (flag_log = 1)
    (@isdefined flag_figure) || (flag_figure = 1)
    (@isdefined flag_outfile) || (flag_outfile = 1)

    # Restrict Weep, REMOVE BEFORE MERGE !!!
    (@isdefined Nrestrict) || (Nrestrict = nothing)
    (@isdefined wndRestrict) || (wndRestrict = nothing)


    ### Check if input files exist ###
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


    ### Write input to struct/dictionary ###
    # rewrite s.t. @isdefined is checked in struct?
    # Use substructures, e.g. struct a2f, ...
    # remove types in definition ::Type
    


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

end

"""
    checkInputFiles()

Check which input files (a2f, dos, weep) exist. Overwrite flags if neccessary.

# Examples
"""
function checkInputFiles(cDOS_flag, include_Weep)
    if ~@isdefined(a2f_file)
        error("a2f-file not specified!")
        
    elseif ~@isdefined(dos_file) && ~@isdefined(dosW_file) && (cDOS_flag == 0 || include_Weep == 1)
        cDOS_flag = 1 
        include_Weep = 0
        print(@yellow "WARNING: ")
        print("No Dos-file specified! Calculating Tc within constant Dos approximation using Anderson-pseudopotential instead\n\n")
    
        if flag_log == 1
            print(log_file, "WARNING: No Dos-file specified! Calculating Tc within constant Dos approximation using Anderson-pseudopotential instead\n\n")
        end
    
    elseif (~(@isdefined(Weep_file) || ~@isdefined(dosW_file))) && include_Weep == 1
        include_Weep = 0
        print(@yellow "WARNING: ")
        print("No Weep/DosWeep-file specified! Calculating Tc using Anderson-pseudopotential instead\n\n")
    
        if flag_log == 1
            print(log_file, "WARNING: No Weep/DosWeep-file specified! Calculating Tc using Anderson-pseudopotential instead\n\n")
        end

    end

    return cDOS_flag, include_Weep

end

"""
    readIn_a2f(a2f_file, indSmear [,unit, nheader, nfooter, nsmear])   

Read in a2f file to solve the isotropic Migdal-Eliashberg equations

The first column must contain the energies. From the second column onwards a2F values for different smearings
"""
function readIn_a2f(a2f_file, indSmear, unit=nothing, nheader=nothing, nfooter=nothing, nsmear=nothing)   
    ### Read in a2f file ###
    a2f_data = readdlm(a2f_file);

    ### Define defaults
    (~isnothing(nheader)) || (nheader = findfirst(isa.(a2f_data[:,1], Number))-1)
    (~isnothing(nfooter)) || (nfooter = size(a2f_data, 1) - findlast(isa.(a2f_data[:,1], Number)))
    (~isnothing(nsmear)) || (nsmear = length(a2f_data[nheader, isa.(a2f_data[nheader,:], Number)])-1)

    ### Remove header & footer
    header = join(a2f_data[1:nheader,:], " ")
    a2f_data = Float64.(a2f_data[nheader+1:end-nfooter, 1:nsmear+1])  # previous version: nheader+1:end-nfooter

    ### Convert omega ###
    omega_raw = a2f_data[:, 1]
    (~isnothing(unit)) || (unit = getUnit(header, "a2F"))
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
function readIn_Dos(dos_file, colFermi=0, spin=1, unit=nothing, nheader=nothing, nfooter=nothing)
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
    (~isnothing(nheader)) || (nheader = findfirst(isa.(dos_data[:,1], Number))-1)
    (~isnothing(nfooter)) || (nfooter = size(dos_data,1) - findlast(isa.(dos_data[:,1], Number)))

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
    (~isnothing(unit)) || (unit = getUnit(header, "Dos"))
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


# Read in Weep
function readIn_Weep(Weep_file, unit=nothing, nheader=nothing, nfooter=nothing)
    """
    Read in Weep file to solve the isotropic Migdal-Eliashberg equations
    All quantities are converted to meV

    -------------------------------------------------------------------
    Input:
        Weep_file:  path to the Weep file
        nheader:    number of header lines in Weep file
        nfooter:    number of footer lines in Weep file
        unit:       units used in Weep file  
                        - "meV"
                        - "eV"
        

    --------------------------------------------------------------------
    Output:
        Weep:       Matrix containing full screened coloumb interaction, 
                    [meV]
    
    --------------------------------------------------------------------
    Comments:
        - Currently Weep data must be in column 3 !!!
    --------------------------------------------------------------------
    """


    ### Read in Weep file ###
    Weep_data = readdlm(Weep_file);

    ### Default values ###
    (~isnothing(nheader)) || (nheader = findfirst(isa.(Weep_data[:,1], Number))-1)
    (~isnothing(nfooter)) || (nfooter = size(Weep_data,1) - findlast(isa.(Weep_data[:,1], Number)))

    ### Remove header & footer
    header = join(Weep_data[1:nheader,:], " ")
    Weep = Float64.(Weep_data[nheader+1:end-nfooter, 3])
    
    ### reshape to matrixs
    Weep = transpose(reshape(Weep, (Int(sqrt(size(Weep, 1))), Int(sqrt(size(Weep, 1))))))

    ### remove outliers from Weep ###
    outlier = findall(Weep .< 0)    #.|| Weep .> 1e5	
    for k in eachindex(outlier)
        Weep[outlier[k]] = 0
    end

    ### Convert ###
    (~isnothing(unit)) || (unit = getUnit(header, "Weep"))
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


# Remove zeros Dos
function neglectZeros(Dos::Vector{Float64}, energies::Vector{Float64}, Weep::Matrix{Float64})
    """
    Remove the leading zeros in the Dos
    There will be no contribution if the dos = 0

    -------------------------------------------------------------------
    Input:
        Dos:        Dos data
        energies:   energy grid points
        Weep:       Weep data

    --------------------------------------------------------------------
    Output:
        -
    
    --------------------------------------------------------------------
    Comments:
        - 
    --------------------------------------------------------------------
    """


    idxLower = findfirst(Dos .!= 0)
    idxUpper = findlast(Dos .!= 0)


    return     Dos[idxLower:idxUpper], energies[idxLower:idxUpper], Weep[idxLower:idxUpper, idxLower:idxUpper]

end

function neglectZeros(Dos::Vector{Float64}, energies::Vector{Float64})
    """
    Remove the leading zeros in the Dos
    There will be no contribution if the dos = 0

    -------------------------------------------------------------------
    Input:
        Dos:        Dos data
        energies:   energy grid points
        Weep:       Weep data

    --------------------------------------------------------------------
    Output:
        -
    
    --------------------------------------------------------------------
    Comments:
        - 
    --------------------------------------------------------------------
    """


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