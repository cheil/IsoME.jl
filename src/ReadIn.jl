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