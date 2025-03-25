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

        ### Imaginary Eliashberg
        if isempty(inp.path_selfEnergy)
            eliash = transfer_fields(inp, arguments)  # Transfer matching fields to arguments

            # set defaults
            eliash.flag_writeSelfEnergy = 1
            eliash.cDOS_flag            = 1
            eliash.include_Weep         = 0

            # solve eliashberg 
            EliashbergSolver(eliash)

            # self energy path
            inp.path_selfEnergy = inp.outdir*"SelfEnergy/"
            inp.muc_ME = eliash.muc_ME
        else

            # read muc_ME, temps,... from Info file in path_selfEnergy?

            inp, log_file, errorLogger = createDirectory(inp, strIsoME)
        end

"""
 Not finishted
"""
        ### read inputs
        matval = ()
        console = Dict()
        try
            inp, console, matval, ML_Tc = InputParser(inp, log_file)        ### !!!! Not finished !!! --> read in a2f and self energy components at each temperature
        catch ex
            # crash file
            writeToCrashFile(inp)

            # console / log file
            printError("while reading the inputs. Stopping now!", ex, log_file, errorLogger)
 
            rethrow(ex)
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
    end


    # material specific values
    matval = (a2f_omega, a2f, dos_en, dos, Weep, dosef, idx_ef, ndos, BCS_gap, idxShiftcut)

    return inp, console, matval, ML_Tc
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

    # check input files / cDOS & Weep
    if ~isfile(inp.a2f_file)
        error("Invalid path to a2f-file!")
    end

    # 
    if any(inp.temps .< 0)
        error("Invalid temperature specification!")
    end

    return inp

end
