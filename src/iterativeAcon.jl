"""
    File containing the iterative analytic continuation of the electron self-energy as in
        (Marsiglio et. al., Iterative analytic continuation of the electron self-energy to the real axis)

"""


"""
    

Iterative analytic continuation of Migdal-Eliashberg theory to the real frequency axis. 
The imaginary results from a cDOS+μ calculation are reuqired as input.
"""
function iterativeAcon(inp::FINDAGOODNAME)

    dt = @elapsed begin

        strIsoME = printIsoME()

        ### Check input
        inp = checkInput(inp)

        #=
        Two modes:
            1) Eliashberg calculations already finished, requires path to info-file and self-energy components
            2) no eliashberg calculations, start EliashbergSolver() first
        =#

        ### Imaginary Eliashberg
        if isempty(inp.path_selfEnergy) || isempty(inp.path_InfoFile)
            eliash = transfer_fields(inp, arguments)  # Transfer matching fields to arguments

            # set defaults
            eliash.flag_writeSelfEnergy = 1
            eliash.cDOS_flag            = 1
            eliash.include_Weep         = 0

            # solve eliashberg 
            EliashbergSolver(eliash)

            # self energy path
            inp.path_selfEnergy = inp.outdir*"SelfEnergy/"
            inp.path_Info = eliash.outdir
            inp.muc_ME = eliash.muc_ME
            inp.temps = eliash.temps
        else
            inp, log_file, errorLogger = createDirectory(inp, strIsoME)

            ### Read Info file ###
            inp.temps, inp.muc_ME = readInfoFile(inp.path_InfoFile)     # is it safe to assume that info is path_selfEnergy/..?
        end

"""
 Not finished
"""
        ### read inputs
        matval = ()
        console = Dict()
        try
            inp, console, matval = InputParser(inp, log_file)        ### !!!! Not finished !!! --> read in a2f and self energy components at each temperature
        catch ex
            # crash file
            writeToCrashFile(inp)

            # console / log file
            printError("while reading the inputs. Stopping now!", ex, log_file, errorLogger)
 
            rethrow(ex)
        end


        ### loop over temperatures
        inp.temps = sort(inp.temps)
        for itemp in inp.temps 
            solveIterativeAcon(itemp, inp, matval)

        end
        

    end

    # print time elapsed
    print("\nTotal Runtime: ", round(dt, digits=2), " seconds\n")
    # log file
    print(log_file, "\nTotal Runtime: ", round(dt, digits=2), " seconds\n")
    
    # close & save
    close(log_file)


end


"""
    InputParser(inp, log_file)

Read, convert, preprocess inputs for iterative_ACON().
"""
function InputParser(inp::FINDAGOODNAME, log_file)

    ### Init table size ###
    console = Dict()

    # header table
    console["header"] = ["it", "znormi", "deltai", "err_delta"]
    #width table
    console["width"] = [8, 14, 14, 14]
    # precision data
    console["precision"] = [0, 4, 4, 5]
    

    console = formatTableHeader(console)

    console = printStartMessage(console, log_file)
  

    ########## READ-IN ##########
    ####### a2f file #######
    a2f_omega, a2f, inp.ind_smear, inp.a2f_unit = readIn_a2f(inp.a2f_file, inp.ind_smear, inp.a2f_unit, inp.nheader_a2f, inp.nfooter_a2f, inp.nsmear)

    # material specific values
    matval = (a2f_omega, a2f)

    return inp, console, matval
end


"""

Read the info.dat file
"""
function readInfoFile(path_InfoFile)


    data = readdlm(path_InfoFile)
    muc_ME = Float64(data[findfirst(data[:,1].=="muc_ME:"),2])
    temps = numbersFromString(join(data[findfirst(data[:,1].=="temps:"),:]))

    return temps, muc_ME
    
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

function readSelfEnergy(itemp, path)
    
    data = readdlm(path*"Delta_"*string(itemp)*"K.dat")
    wsi = Float64.(data[2:end,1])
    Deltai = Float64.(data[2:end,2])

    data = readdlm(path*"Z_"*string(itemp)*"K.dat")
    Zi = Float64.(data[2:end,2])

    return wsi, Zi, Deltai

end

"""

Solve the iteravie ACON at the given temperature
"""
function solveIterativeAcon(itemp, inp, matval)
    # destruct inputs
    (a2f_omega, a2f) = matval
    (; outdir) = inp

    ### Read imaginary self-energy ###
    wsi, Zi, Deltai = readSelfEnergy(itemp, inp.path_selfEnergy)

    ### real axis ###
    ws = range(-inp.real_c, inp.real_c, inp.numReal_c)

    # lambda, it is faster to calcualte lambda each time
    lambda = calcLambda(ws, wsi, a2f, a2f_omega)            # !!! VERY SLOW !!!

    ### cDOS + mu case ###
    iterativeAconEq(itemp, ws, wsi, Zi, Deltai, λ)


end

"""

Eq. (3.5a), (3.5b) in "Iterative analytic continuation of the electron self-energy to the real energy axis", Marsiglio et.al, (1987)
"""
function iterativeAconEq(itemp, ws, wsi, Zi, Deltai, λ)

    ZimagSum = ws + i*π*kb*itemp*sum(ws./sqrt.(ws.^2 .- (Deltai^2)') .* () )


end


"""

calculate lambda with complex argument
lambda[ws, wsi]
"""
function calcLambda(ws::StepRangeLen{Float64}, wsi::Vector{Float64}, a2f::Vector{Float64}, a2f_omega::StepRangeLen{Float64})

    lambda = Array{ComplexF64}(undef, length(ws), length(wsi))
    for (omega, idx) in zip(ws, collect(1:length(ws)))
        z = omega .+ wsi
        integrand = 2 * a2f' .* a2f_omega' ./ (z .^ 2 .+ a2f_omega' .^ 2)
        lambda[idx, :] = trapz(a2f_omega, integrand)
    end

    return lambda
end


""" 

Transfer commmon fields from one struct to another
"""
function transfer_fields(from_struct, ::Type{T}) where T

    struct2tuple(p) = (; (v=>getfield(p, v) for v in fieldnames(typeof(p)))...)

    from_nt = struct2tuple(from_struct)  # Convert struct to named tuple
    common_fields = filter(kv -> hasfield(T, kv[1]), pairs(from_nt))  # Keep matching fields
    return T(; common_fields...)  # Pass them to construct T
end


"""
    checkInput!(inp)

Check if a2F-file exists and temps are specified.
"""
function checkInput(inp::FINDAGOODNAME)

    # check input files
    if ~isfile(inp.a2f_file)
        error("Invalid path to a2f-file!")
    end

    # 
    #if any(inp.temps .< 0)
    #    error("Invalid temperature specification!")
    #end

    return inp

end
