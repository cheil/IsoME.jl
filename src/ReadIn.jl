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
    InputParser!(inp, log_file)

Read, convert, preprocess inputs for EliashbergSolver() (mode=0) and RealAxisSolver() (mode=1).
"""
function InputParser!(inp::arguments, log_file; mode::Int64=0)
    # mode = 0: EliashbergSolver()
    # mode = 1: RealAxisSolver()

    ### Init table size ###
    console = initOutputTable(inp, mode=mode)

    console = formatTableHeader(console)

    console = printStartMessage(console, inp, log_file, mode=mode)


    # ================ READ-IN ================ #
    # --------------- a2f file ---------------- #
    a2f_omega, a2f, a2f_itp, inp.ind_smear, inp.a2f_unit = readIn_a2f(inp.a2f_file, inp.ind_smear, inp.a2f_unit, inp.nheader_a2f, inp.nfooter_a2f, inp.nsmear)


    # ------------- Dos and Weep -------------- #
    if isfile(inp.dos_file) && (inp.cDOS_flag == 0 || (isnan(inp.mu) && isfile(inp.Weep_file)))

        # read dos
        dos_en, dos, ef, inp.dos_unit = readIn_Dos(inp.dos_file, inp.ef, inp.spinDos, inp.dos_unit, inp.nheader_dos, inp.nfooter_dos, outdir=inp.outdir, logFile=log_file)
        inp.ef = ef

        # remove zeros at begining/end of dos
        dos, dos_en = discardZeros(dos, dos_en)

        ### Interpolation ###
        # Interpolation Object DoS
        epsilonItp = range(dos_en[1], dos_en[end], length(dos_en))
        itpDos = scale(interpolate(dos, BSpline(Linear())), epsilonItp)

        # encut in meV
        if (-inp.encut < dos_en[1]) || (inp.encut > dos_en[end])
            text = "Energy cutoff exceeds range of DOS!"
            printWarning(text, log_file)
        end

        if inp.include_Weep == 1 || ((isfile(inp.Weep_file) && (isempty(inp.Wen_file) || isfile(inp.Wen_file))) && isnan(inp.mu))
            # read Weep + energy grid points
            Weep, Wen, inp.efW, inp.Weep_unit = readIn_Weep(inp.Weep_file, inp.Wen_file, inp.Weep_col, inp.Wen_col, inp.efW, inp.Weep_unit, inp.nheader_Weep, inp.nfooter_Weep, inp.nheader_Wen, inp.nfooter_Wen, outdir=inp.outdir, logFile=log_file)

            ### Interpolation ###
            # Interpolation Object Weep
            epsilonItp = range(Wen[1], Wen[end], length(Wen))       
            itpWeep = scale(interpolate(Weep, BSpline(Linear())), (epsilonItp, epsilonItp)) 

            if mode == 0 
                # interpolate 
                dos_en, dos, Weep = interpolateInputs(itpDos, dos_en, inp.itpStepSize, inp.itpBounds, inp.encut, itpWeep=itpWeep, Wen=Wen)
            else
                # interpolate dos/Weep on same grid
                lowCut  = max(dos_en[1],   -inp.encut)
                highCut = min(dos_en[end],  inp.encut)
                enStart, enEnd = gridWindow(Wen, lowCut, highCut, "W-energy")
                de = inp.depsilon
                dos_en = collect(enStart:de:enEnd)
                dos = itpDos(dos_en)
                Weep = itpWeep(dos_en, dos_en)

            end

            # mu, idxWef = idx_ef
            (!isnan(inp.mu)) || (idxWef = findmin(abs.(dos_en))[2]; inp.mu = dos[idxWef] .* Weep[idxWef, idxWef])

        else
            if mode == 0
                # interpolate 
                dos_en, dos, Weep = interpolateInputs(itpDos, dos_en, inp.itpStepSize, inp.itpBounds, inp.encut)
            else 
                enStart, enEnd = gridWindow(dos_en, -inp.encut, inp.encut, "DOS-energy")
                de = inp.depsilon
                dos_en = collect(enStart:de:enEnd)
                dos = itpDos(dos_en)
                Weep = Matrix{Float64}(undef, 0, 0)
            end
        end

    else
        # default values (typed empties so matval stays concrete: Vector{Float64}, not Vector{Any})
        dos = Float64[]
        dos_en = Float64[]
        ef = NaN
        Weep = Matrix{Float64}(undef, 0, 0)
    end

    if inp.include_Weep == 0 && inp.cDOS_flag == 1
        # default values in case no dos-file is given
        ndos = -1
        idx_ef = -1
        dosef = -1.0        # Float64 to match dos[idx_ef] in the else-branch (keeps matval[6] concrete)
    else
        # length energy vector
        ndos = size(dos_en, 1)

        # index of fermi energy
        idx_ef = findmin(abs.(dos_en))
        idx_ef = idx_ef[2]

        # dos at ef
        dosef = dos[idx_ef]
    end

    ### optionally rescale α²F to a better N(ε_F) (once, before any consumer: μ*, AD, matval)
    a2f = rescale_a2F(a2f, inp, dosef, log_file)

    ### calc mu*'s
    # both solvers use the same frequency cutoff
    phonon_cutoff = inp.omega_c

    if isnan(inp.mu) && isnan(inp.muc_ME) && isnan(inp.muc_AD)
        # defaut value
        inp.muc_AD = 0.12
        calcMucME(inp, a2f, a2f_omega, phonon_cutoff, log_file)

    elseif isnan(inp.muc_AD) && isnan(inp.muc_ME)
        if !isnan(inp.typEl)
            calcMucs(inp, inp.typEl, a2f, a2f_omega, phonon_cutoff, log_file)

        elseif ~(isnan(ef) || ef == 0)
            inp.typEl = ef
            calcMucs(inp, inp.typEl, a2f, a2f_omega, phonon_cutoff, log_file)

        elseif ~(isnan(inp.efW) || inp.efW == 0)
            inp.typEl = inp.efW
            calcMucs(inp, inp.typEl, a2f, a2f_omega, phonon_cutoff, log_file)

        else
            text = "Unable to calculate μ* from μ without a typical electron energy!"
            text *= "\nConsider setting typEl, ef or efW manually."
            text *= "\nUsing μ*_AD = 0.12 instead."
            text *= "\nSee the pseudopotential section of the Input documentation for the μ → μ* relation."
            printWarning(text, log_file)

            inp.muc_AD = 0.12
            calcMucME(inp, a2f, a2f_omega, phonon_cutoff, log_file)
        end

    elseif inp.include_Weep == 0 && !isnan(inp.muc_AD) && isnan(inp.muc_ME)
        calcMucME(inp, a2f, a2f_omega, phonon_cutoff, log_file)

    elseif !isnan(inp.muc_ME) && isnan(inp.muc_AD)
        calcMucAD(inp, a2f, a2f_omega, phonon_cutoff)
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
    matval = (a2f_omega, a2f, dos_en, dos, Weep, dosef, idx_ef, ndos, BCS_gap)

    return console, matval, ML_Tc
end



"""
    initOutputTable(inp; mode=0)

initialize output table: header names, column width and precision
"""
function initOutputTable(inp::arguments; mode::Int64=0)#
    # mode = 0: EliashbergSolver()
    # mode = 1: RealAxisSolver()

    # Both modes access the table via console.cDOS / console.vDOS. The imaginary-axis
    # solver (mode 0) populates only the slot selected by cDOS_flag; the other stays empty.
    console = Console()
    if mode == 0
        if inp.include_Weep == 1 && inp.cDOS_flag == 0
            console.vDOS = TableSpec(["it", "phic", "phiph", "znormi", "shifti", "ef-mu", "deltai", "err_delta"], [8, 10, 10, 10, 10, 10, 10, 11], [0, 2, 2, 2, 2, 2, 2, 5])

        elseif inp.include_Weep == 1 && inp.cDOS_flag == 1
            console.cDOS = TableSpec(["it", "phic", "phiph", "znormi", "deltai", "err_delta"], [8, 10, 10, 10, 10, 11], [0, 2, 2, 2, 2, 5])

        elseif inp.include_Weep == 0 && inp.cDOS_flag == 0
            console.vDOS = TableSpec(["it", "znormi", "shifti", "ef-mu", "deltai", "err_delta"], [8, 10, 10, 11, 10, 11], [0, 2, 2, 2, 2, 5])

        elseif inp.include_Weep == 0 && inp.cDOS_flag == 1
            console.cDOS = TableSpec(["it", "znormi", "deltai", "err_delta"], [8, 14, 14, 14], [0, 4, 4, 5])

        else
            error("Unknown mode! Check if the cDOS_flag and include_Weep flag are set correctly!")
        end

    elseif mode == 1

        if inp.include_Weep in (0, 1)
            console.cDOS = TableSpec(["it", "Re(Z)", "Im(Z)", "Re(Δ)", "Im(Δ)", "Error Δ"], [8, 10, 10, 10, 10, 11], [0, 4, 4, 4, 4, 5])
            console.vDOS = TableSpec(["it", "Re(Z)", "Im(Z)", "Re(χ)", "Im(χ)", "ef-mu", "Re(Δ)", "Im(Δ)", "Error Δ"], [8, 10, 10, 10, 10, 10, 10, 10, 11], [0, 4, 4, 4, 4, 2, 4, 4, 5])
        end

    end

    return console
end

"""
    createDirectory!(inp, strIsoME)

Create directory.
"""
function createDirectory!(inp::arguments, strIsoME::String)

    # the logger in place before the run, handed back so that the solver can put it back
    # when it is done (see `finalizeLogging`): IsoME must not leave a session logging into
    # a closed log file of a finished run
    prevLogger = global_logger()

    if inp.testMode # hidden outputs in test mode
        log_file = IOBuffer()
        errorLogger = SimpleLogger(log_file, Logging.Error)
    else
        # create output directory
        if isempty(inp.outdir)
            inp.outdir = "./"
        elseif ~(inp.outdir[end] == '/' || inp.outdir[end] == '\\')
            inp.outdir = inp.outdir * "/"
        end

        # never write into an existing directory: append a run counter instead, so that
        # consecutive runs into the same path land in <outdir>_1, <outdir>_2, ...
        idxDir = 1
        tempDir = inp.outdir[1:end-1]
        while isdir(inp.outdir)
            inp.outdir = tempDir * "_" * string(idxDir) * "/"
            idxDir += 1
        end

        try
            mkpath(inp.outdir)
        catch ex
            ex isa InterruptException && rethrow(ex)
            error("Couldn't write into " * inp.outdir * "! Outdir may not writable or an invalid path.\n\n")
        end

        log_file = open(inp.outdir * "log.txt", "w")
        print(log_file, strIsoME)

        # the log file is buffered: make sure what has been written so far reaches
        # the disk even if the process is ended without unwinding the solver.
        # capture an alias that is assigned once, so that `log_file` - which is set in
        # both branches above - is not boxed and keeps its type in the return value
        let lf = log_file
            atexit(() -> flushLog(lf))
        end

        # logging to console and log-file (@warn,...)
        errorLogger = SimpleLogger(log_file, Logging.Error)
        file_logger = SimpleLogger(log_file)
        tee_logger = TeeLogger(ConsoleLogger(), file_logger)
        global_logger(tee_logger)
    end

    return log_file, errorLogger, prevLogger

end


"""
    finalizeLogging(log_file, prevLogger)

End of a run: flush and close the log file and put the logger that was active before the
run back in place. Never throws.

`createDirectory!` redirects `@warn`/`@info` into the run's log file through the global
logger. Without this, that redirection would outlive the run and every later message of the
session would be written to a log file that has since been closed.
"""
function finalizeLogging(log_file, prevLogger)
    closeLog(log_file)
    try
        global_logger(prevLogger)
    catch
    end
    return nothing
end


"""
    checkInput!(inp; realSolver=false)

Check which input files (a2f, dos, weep) exist, that the mode flags are valid, and that
`encut` stays inside `omega_c`.
For real axis solver: Additionally check if cDOS+W mode has been chosen
"""
function checkInput!(inp::arguments; realSolver::Bool=false)

    # re-run the range check on the struct: the constructor already rejected a negative
    # value, but the fields are mutable and may have been assigned since
    _check_negative(k => getfield(inp, k) for k in keys(_nonNegative))

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

    if inp.cDOS_flag ∉ (0,1)
        error("Invalid cDOS_flag value. Use 0 for a variable density of states (vDOS) or 1 for a constant density of states (cDOS).")
    end

    if inp.include_Weep ∉ (0, 1) 
        error("Invalid include_Weep value. Use 0 for the μ approximation or 1 for the W(ε, ε′) interaction.")
    end

    if realSolver && inp.include_Weep == 1 && inp.cDOS_flag == 1
        error("The real-axis solver only supports the vDOS+W approximation (include_Weep = 1 requires cDOS_flag = 0). Set cDOS_flag = 0.\n\n")
    end

    # vDOS, both axes: the frequency sum/integral of the μ-update and of χ only shows its
    # correct limiting behaviour as long as the frequency reaches at least as far as ε, so
    # the ε-window must stay inside the frequency cutoff. Require encut ≤ omega_c and clamp
    # encut down. On the imaginary axis omega_c is the Matsubara cutoff, on the real axis
    # the cutoff of the ω-grid; the requirement is the same in both cases.
    if inp.cDOS_flag == 0 && inp.encut > inp.omega_c
        @warn "encut = $(inp.encut) exceeds omega_c = $(inp.omega_c); the ε-integration has to stay inside the frequency cutoff. Setting encut = $(inp.omega_c). Increase omega_c to keep the larger ε-window. See the μ-update section of the Troubleshooting page."
        inp.encut = inp.omega_c
    end

    return nothing

end

"""
    rescale_a2F(a2f, inp, dosef, log_file) -> a2f

Optionally rescale α²F to a better density of states at the Fermi level.

If `inp.a2f_Nef` is set (not `NaN`), α²F is multiplied by `a2f_Nef / N(ε_F)`, where `a2f_Nef` is the
N(ε_F) that the α²F calculation used and `N(ε_F)` is taken from the DOS file, so that α²F and the
electronic normalization used in the Eliashberg equations share the *same* N(ε_F). 

The better N(ε_F) is taken from `dosef` when the read-in already produced it (vDOS / Weep+μ). In
cDOS+μ, where no DOS is read, it is read from `dos_file` if one is given (a local read for the
rescale only; the cDOS `matval` sentinels are untouched). Without any DOS file the rescale is
skipped with a warning. No-op when `a2f_Nef` is `NaN`.
"""
function rescale_a2F(a2f, inp, dosef, log_file)
    isnan(inp.a2f_Nef) && return a2f                        # feature off

    inp.a2f_Nef > 0 || error("a2f_Nef must be positive (got $(inp.a2f_Nef)).")

    # better N(ε_F): use the read-in dosef if available, otherwise read it from the DOS file (cDOS+μ)
    Nef = dosef
    if !(Nef > 0) && isfile(inp.dos_file)
        dos_en_r, dos_r, _, _ = readIn_Dos(inp.dos_file, inp.ef, inp.spinDos, inp.dos_unit, inp.nheader_dos, inp.nfooter_dos, outdir=inp.outdir, logFile=log_file)
        dos_r, dos_en_r = discardZeros(dos_r, dos_en_r)
        Nef = dos_r[findmin(abs.(dos_en_r))[2]]             # N(ε_F): DOS at the energy closest to ε_F
    end

    if !(Nef > 0)
        printWarning("a2f_Nef is set but no DOS file is available: α²F cannot be rescaled without a DOS file. Ignoring a2f_Nef.", log_file)
        return a2f
    end

    factor = inp.a2f_Nef / Nef
    printTee(log_file, "\nα²F rescaled to N(ε_F) from the DOS file: a2f_Nef = " * string(round(inp.a2f_Nef, sigdigits=5)) *
                       ", N_dos(ε_F) = " * string(round(Nef, sigdigits=5)) *
                       ", factor = " * string(round(factor, sigdigits=5)) * "\n")
    return a2f .* factor
end


"""
    autoHeaderFooter(data, nheader, nfooter, nameFile) -> nheader, nfooter

Fill in the header/footer size of a read-in file where it was left to auto-detection
(-1): everything above the first numeric entry of the first column is header, everything
below the last one is footer.
"""
function autoHeaderFooter(data::Matrix{Any}, nheader::Int, nfooter::Int, nameFile::AbstractString)
    (nheader >= 0 && nfooter >= 0) && return nheader, nfooter

    numeric = isa.(data[:, 1], Number)
    firstNum = findfirst(numeric)
    isnothing(firstNum) && error("The first column of the " * nameFile * "-file holds no numeric entry, so header and footer can not be detected. Check the file and its column layout, or set the header/footer size manually.\n\n")

    (nheader >= 0) || (nheader = firstNum - 1)
    # a first numeric entry guarantees a last one
    (nfooter >= 0) || (nfooter = size(data, 1) - findlast(numeric)::Int)

    return nheader, nfooter
end


"""
    readIn_a2f(a2f_file, indSmear=-1, unit="", nheader=-1, nfooter=-1, nsmear=-1)

Read in a2f file to solve the isotropic Migdal-Eliashberg equations

The first column must contain the energies, the second column onwards a2F values for different smearings
"""
function readIn_a2f(a2f_file, indSmear::Int=-1, unit="", nheader::Int=-1, nfooter::Int=-1, nsmear::Int=-1)
    ### Read in a2f file ###
    # `Any` forces a Matrix{Any}: without it readdlm returns a Matrix{Float64} for a file
    # without a header and a Matrix{Any} for one with, and the resulting union makes the
    # whole header/footer detection below type-unstable
    a2f_data = readdlm(a2f_file, Any)::Matrix{Any}

    ### Define defaults (-1 == auto-detect)
    nheader, nfooter = autoHeaderFooter(a2f_data, nheader, nfooter, "a2F")
    # ::Int assert: the raw (Any) a2f_data makes length() infer Any, which would leave
    # nsmear/indSmear (and hence the returned a2f) type-unstable.
    (nsmear >= 0)   || (nsmear = (length(a2f_data[nheader+1, isa.(a2f_data[nheader+1, :], Number)]) - 1)::Int)
    (indSmear >= 0) || (indSmear = Int64(ceil(nsmear / 2)))

    ### Remove header & footer
    header = join(a2f_data[1:nheader, :], " ")
    # type-assert to a concrete Matrix{Float64}: readdlm infers as Any, so without this
    # the whole numeric pipeline (and the return) stays Any-typed.
    a2f_data = Float64.(a2f_data[nheader+1:end-nfooter, 1:nsmear+1])::Matrix{Float64}

    ### Convert omega ###
    omega_raw = a2f_data[:, 1]
    (~isempty(unit)) || (unit = getUnit(header, "a2F"))
    if "meV" == unit
        omega_raw = omega_raw
    elseif "eV" == unit
        omega_raw = omega_raw * 1000
    elseif "THz" == unit
        omega_raw = omega_raw * THz2meV
    elseif "Ry" == unit
        omega_raw = omega_raw * Ry2meV
    elseif "Ha" == unit     # Hartree
        omega_raw = omega_raw .* Ry2meV * 2
    else
        error("Could not determine the unit of the a2F-file. Set it manually via a2f_unit or check the file header (supported units: meV, eV, THz, Ry, Ha).")
    end

    ### a2f for one smearing ###
    a2f_raw = a2f_data[:, indSmear+1]

    ### interpolate a2F on 10x finer grid ###
    omega = range(1e-2, stop=omega_raw[end], length=size(omega_raw)[1])   
    a2f_itp = linear_interpolation(omega_raw, a2f_raw, extrapolation_bc=0)     
    a2f = a2f_itp(omega)
    a2f[a2f.<0.0] .= 0.0
    
    return omega, a2f, a2f_itp, indSmear, unit

end


# Read in DOS
"""
    readIn_Dos(dos_file, ef=NaN, spin=2, unit="", nheader=-1, nfooter=-1; outdir ="./", logFile = nothing)

Read the dos file. All quantities are converted to meV.

The energies must be in column 1 and the dos in column 2
"""
function readIn_Dos(dos_file, ef::Float64=NaN, spin=2, unit="", nheader::Int=-1, nfooter::Int=-1; outdir="./", logFile=nothing)

    ### Read in dos file ###
    dos_data = readdlm(dos_file, Any)::Matrix{Any}

    ### Default values (-1 == auto-detect) ###
    nheader, nfooter = autoHeaderFooter(dos_data, nheader, nfooter, "Dos")

    ### Remove header & footer
    header = dos_data[1:nheader, :]
    dos = Float64.(dos_data[nheader+1:end-nfooter, 1:2])::Matrix{Float64}
    (~isempty(unit)) || (unit = getUnit(join(header, " "), "Dos"))

    ### extract energies and dos
    energies = dos[:, 1]
    dos = dos[:, 2]
    dos[dos.<0.0] .= 0.0 # set negative dos to 0

    ### spin
    dos = dos / spin

    ### Fermi energy
    if isnan(ef)
        ef = extractFermiEnergy(header, unit, "Dos", outdir=outdir, logFile=logFile)
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
    elseif "Ha" == unit     # Hartree
        # the DOS is per unit energy, so it is divided by the *same* factor the energies
        # are multiplied with: 1 Ha = 2 Ry
        energies = energies .* (Ry2meV * 2)
        dos = dos ./ (Ry2meV * 2)
    else
        error("Could not determine the unit of the DOS-file. Set it manually via dos_unit or check the file header (supported units: meV, eV, THz, Ry, Ha).")
    end

    ### Shift energies by ef for cDos ###
    energies = energies .- ef

    return energies, dos, ef, unit
end


"""
    readIn_Weep(Weep_file, Wen_file="", Weep_col=3, Wen_col=1, ef=NaN, unit = "", nheader=-1, nfooter=-1,  nheaderWen=-1, nfooterWen=-1; outdir = "./", logFile = nothing)

Read in Weep file containing the sreened coulomb interaction.
Weep data must be in column 3
"""
function readIn_Weep(Weep_file, Wen_file="", Weep_col=3, Wen_col=1, ef::Float64=NaN, unit="", nheader::Int=-1, nfooter::Int=-1, nheaderWen::Int=-1, nfooterWen::Int=-1; outdir="./", logFile=nothing)

    ### Read in Weep file ###
    Weep_data = readdlm(Weep_file, Any)::Matrix{Any}

    # Default values (-1 == auto-detect)
    nheader, nfooter = autoHeaderFooter(Weep_data, nheader, nfooter, "Weep")

    # Remove header & footer
    header = Weep_data[1:nheader, :]
    Weep = Float64.(Weep_data[nheader+1:end-nfooter, Weep_col])::Vector{Float64}

    # reshape to matrix (Matrix, not a lazy Transpose, so matval carries a concrete Matrix{Float64})
    numWens = Int(sqrt(size(Weep, 1)))
    Weep = Matrix(transpose(reshape(Weep, numWens, numWens)))

    # remove outliers from Weep 
    Weep[Weep.<0] .= 0

    # unit
    (~isempty(unit)) || (unit = getUnit(join(header, " "), "Weep"))

    ### Fermi energy
    if isnan(ef)
        ef = extractFermiEnergy(header, unit, "Weep", outdir=outdir, logFile=logFile)
    end

    ### read in W energies ###
    if isempty(Wen_file)
        #Wen = Float64.(Weep_data[nheader+1:numWens+nheader, Wen_col])
        Wen = unique(Float64.(Weep_data[nheader+1:end-nfooter, Wen_col])::Vector{Float64})
    else
        Wen = readIn_Wen(Wen_file, Wen_col, nheaderWen, nfooterWen)
    end

    ### Convert ###
    if "meV" == unit    # meV
        Weep = Weep
        Wen = Wen
    elseif "eV" == unit     # eV
        Weep = Weep .* 1000
        Wen = Wen .* 1000
    elseif "Ry" == unit      # Ry
        Weep = Weep .* Ry2meV
        Wen = Wen .* Ry2meV
    elseif "Ha" == unit     # Hartree
        Weep = Weep .* Ry2meV * 2
        Wen = Wen .* Ry2meV * 2
    elseif unit == "THz"
        Weep = Weep .* THz2meV
        Wen = Wen .* THz2meV
    else
        error("Could not determine the unit of the Weep-file. Set it manually via Weep_unit or check the file header (supported units: meV, eV, THz, Ry, Ha).")
    end

    # shift by ef
    Wen = Wen .- ef

    return Weep, Wen, ef, unit

end

"""
    readIn_Wen(Wen_file, Wen_col, nheader=nothing, nfooter=nothing)

Read in Wen
Energy grid for Weep
"""
function readIn_Wen(Wen_file, Wen_col, nheader::Int=-1, nfooter::Int=-1)
    ### Read in Weep file ###
    Wen_data = readdlm(Wen_file, Any)::Matrix{Any}

    ### Default values (-1 == auto-detect) ###
    nheader, nfooter = autoHeaderFooter(Wen_data, nheader, nfooter, "Wen")

    ### Remove header & footer
    Wen = Float64.(Wen_data[nheader+1:end-nfooter, Wen_col])::Vector{Float64}

    return Wen

end


# Fermi energy
"""
    extractFermiEnergy(header, unit, nameFile=nothing)

Extract the fermi energy from the header of the input files
"""
function extractFermiEnergy(header, unit, nameFile=nothing; outdir="./", logFile=nothing)

    ef = NaN
    try
        # every branch asserts ::Float64: the header cells are Any, and without it `ef`
        # stays untyped and drags the unit conversion below down with it
        logNums = isa.(header, Number)
        # the cells are Any, so `header .== "="` is not inferred as a Bool mask and drags
        # the sum and every indexing operation below it to Any. Broadcasting a predicate
        # that provably returns Bool keeps the mask - and everything indexed by it - a
        # BitMatrix. It also survives a single-column header, where the mask is empty and
        # `sum` of an empty Matrix{Any} would throw before the keyword search is reached.
        isEqualSign(x)::Bool = x isa AbstractString && x == "="
        maskEq = @view(logNums[:, 2:end]) .& isEqualSign.(@view header[:, 1:end-1])
        if sum(logNums) == 1
            ef = Float64(only(header[logNums])::Number)::Float64
        elseif sum(maskEq) == 1
            ef = Float64(only(@view(header[:, 2:end])[maskEq])::Number)::Float64
        else
            nameFermi = ["efermi", "ef", "fermi", "e_fermi"]

            for name in nameFermi
                # first cell holding the keyword, then the first number at or after it
                idx = findfirst(x -> isa(x, AbstractString) && lowercase(x) == name, header)
                isnothing(idx) && continue

                row, col = Tuple(idx)
                idxNum = findfirst(x -> isa(x, Number), @view header[row, col:end])
                isnothing(idxNum) && error("Found '" * name * "' in the header of the " * nameFile * "-file but no number after it.")

                ef = Float64(header[row, col+idxNum-1]::Number)::Float64
                break
            end
        end

        ### Convert ###
        if "meV" == unit # meV
            ef = ef
        elseif "eV" == unit     # eV
            ef = ef .* 1000
        elseif "THz" == unit    # THz
            ef = ef .* THz2meV
        elseif "Ry" == unit      # Ry
            ef = ef .* Ry2meV
        elseif "Ha" == unit     # Hartree
            ef = ef .* (Ry2meV * 2)
        else
            error("Could not determine the unit of the " * nameFile * "-file. Set it manually via " * nameFile * "_unit or check the file header (supported units: meV, eV, THz, Ry, Ha).")
        end

    catch ex
        ex isa InterruptException && rethrow(ex)
        text = "Error while reading the fermi energy from the " * nameFile * "-file."
        text *= "\nConsider setting the fermi-energy manually (ef or efW) or check the header of the " * nameFile * "-file\n\n"
        error(text)
    end

    # extraction ran without throwing but matched nothing -> ef still NaN
    isnan(ef) && error("Could not extract the Fermi energy from the " * nameFile * "-file.\nSet it manually (ef or efW) or check the file header.\n\n")

    return ef::Float64
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

    # No prompt here: the solver is routinely run non-interactively (batch queue, CI,
    # notebook), where reading from stdin blocks forever or consumes unrelated input.
    field = nameFile == "a2F"  ? "a2f_unit"  :
            nameFile == "Dos"  ? "dos_unit"  :
            nameFile == "Weep" ? "Weep_unit" : "the corresponding *_unit input"
    error("Could not determine the unit of the " * string(nameFile) * "-file from its header.\n" *
          "Set it manually via " * field * " (supported units: meV, eV, THz, Ry, Ha), " *
          "or add the unit to the file header.\n\n")
end


"""
    gridWindow(en, lowCut, highCut, nameGrid) -> enStart, enEnd

First point of the grid `en` strictly above `lowCut` and last one strictly below `highCut`.

Guarded: where the window and the grid do not overlap, `findfirst`/`findlast` return
`nothing` and indexing with it would only surface as a `MethodError` further down. Here it
gives a message naming both ranges instead.
"""
function gridWindow(en::AbstractVector, lowCut::Real, highCut::Real, nameGrid::AbstractString)
    iLow  = findfirst(>(lowCut), en)
    iHigh = findlast(<(highCut), en)

    (isnothing(iLow) || isnothing(iHigh) || iLow > iHigh) && error(
        "The " * nameGrid * " grid (" * string(round(first(en), sigdigits=5)) * " … " *
        string(round(last(en), sigdigits=5)) * " meV) does not overlap the requested energy window [" *
        string(round(lowCut, sigdigits=5)) * ", " * string(round(highCut, sigdigits=5)) *
        "] meV. Check encut and the energy range of the input file.\n\n")

    return en[iLow], en[iHigh]
end


"""
    a2fSupportMax(a2f_omega, a2f) -> ω_max

Largest frequency at which α²F rises above the 1e-2 support threshold, i.e. the
characteristic phonon cutoff entering the μ* conversion formulas.

Guarded: an α²F that stays below the threshold everywhere leaves the masked vector empty,
and `maximum` of an empty collection would abort with an `ArgumentError` that names neither
the file nor the smearing.
"""
function a2fSupportMax(a2f_omega, a2f)
    ω = a2f_omega[a2f.>0.01]
    isempty(ω) && error("α²F stays below 1e-2 over the whole frequency range, so no characteristic phonon frequency can be determined. Check the a2F-file and the selected smearing (ind_smear), or set μ* manually via muc_AD / muc_ME.\n\n")
    return maximum(ω)
end


"""
    discardZeros(Dos, energies)

Discard zeros in dos.
"""
function discardZeros(Dos::Vector{Float64}, energies::Vector{Float64})
    idxLower = findfirst(!iszero, Dos)
    idxUpper = findlast(!iszero, Dos)

    # without this the `nothing` would only surface as a MethodError in the range below
    isnothing(idxLower) && error("The density of states is zero over the whole energy range. Check the dos-file and the column layout (dos_file, spinDos).\n\n")
    isnothing(idxUpper) && error("The density of states is zero over the whole energy range. Check the dos-file and the column layout (dos_file, spinDos).\n\n")

    return Dos[idxLower:idxUpper], energies[idxLower:idxUpper]
end


"""
    calcMucME(inp, console, a2f_omega, log_file)

Calculate μ*_ME from μ*_AD using formula as in 
Pellegrini, Ab initio methods for superconductivity 
DOI: 10.1038/s42254-024-00738-9
"""
function calcMucME(inp, a2f, a2f_omega, phonon_cutoff, log_file)
    inp.muc_ME = inp.muc_AD / (1 + inp.muc_AD * log(a2fSupportMax(a2f_omega, a2f) / phonon_cutoff))

    if inp.muc_ME < 0 || inp.muc_ME > 0.8 || inp.muc_ME > 3 * inp.muc_AD
        inp.muc_ME = minimum([3 * inp.muc_AD, 0.8])

        text = "Couldn't calculate a reasonable μ*_ME from μ*_AD."
        text *= "\nUsing μ*_ME = minimum(3*μ*_AD, 0.8) instead."
        text *= "\nCheck muc_ME.png and consider setting μ* manually or changing the Matsubara cutoff!"
        text *= "\nSee the μ* conversion section of the Troubleshooting page and the pseudopotential section of the Input documentation."
        printWarning(text, log_file)

        wc_plot = range(min(100, floor(phonon_cutoff/2)), max(ceil(2*phonon_cutoff), 1e4), 500)
        muc_ME_plot = inp.muc_AD ./ (1 .+ inp.muc_AD .* log.(a2fSupportMax(a2f_omega, a2f) ./ wc_plot))

        plot(wc_plot, muc_ME_plot)
        savefig(inp.outdir*"muc_ME.png")
    end
end


"""
    calcMucAD(inp, console, a2f_omega)

Calculate μ*_AD from μ*_ME using formula (30) in
Pellegrini, Ab initio methods for superconductivity 
DOI: 10.1038/s42254-024-00738-9
"""
function calcMucAD(inp, a2f, a2f_omega, phonon_cutoff)
    inp.muc_AD = inp.muc_ME / (1 - inp.muc_ME * log(a2fSupportMax(a2f_omega, a2f) / phonon_cutoff))

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
function calcMucs(inp, ef, a2f, a2f_omega, phonon_cutoff, log_file)
    inp.muc_AD = inp.mu / (1 + inp.mu * log(ef / a2fSupportMax(a2f_omega, a2f)))

    # μ*_ME < 4*μ
    if phonon_cutoff > ef * exp(3 / (4 * inp.mu))
        phonon_cutoff = ef * exp(3 / (4 * inp.mu))

        text = "Matsubara cutoff would lead to μ*_ME > 4*μ."
        text *= "\nA smaller cutoff has been used for the μ → μ*_ME conversion; omega_c itself is unchanged,"
        text *= "\nso the solver still runs at omega_c = " * string(inp.omega_c) * " meV."
        text *= "\nCheck muc_ME.png and the typical electronic energy typEl!"
        text *= "\nSee the μ* conversion section of the Troubleshooting page and the pseudopotential section of the Input documentation."
        printWarning(text, log_file)

        wc_plot = range(min(100, floor(phonon_cutoff/2)), max(ceil(2*phonon_cutoff), 1e4), 500)
        muc_ME_plot = inp.muc_AD ./ (1 .+ inp.muc_AD .* log.(a2fSupportMax(a2f_omega, a2f) ./ wc_plot))

        plot(wc_plot, muc_ME_plot)
        savefig(inp.outdir*"muc_ME.png")
    end

    if inp.include_Weep == 0
        inp.muc_ME = inp.mu / (1 + inp.mu * log(ef / phonon_cutoff))
    end

    checkMucsFromMu(inp, ef, a2fSupportMax(a2f_omega, a2f), phonon_cutoff)
end


"""
    checkMucsFromMu(inp, ef, omega_ph, phonon_cutoff)

Reject a negative μ* produced by the μ → μ* conversion in [`calcMucs`](@ref).

Unlike `calcMucME`/`calcMucAD`, which fall back to a sensible value, there is nothing to
fall back to here: a negative μ* means the Morel-Anderson denominator `1 + μ·ln(ε_el/ω)`
has gone through zero, so the conversion is outside its range of validity and any value it
returns is meaningless. It would otherwise enter the equations as an *attractive* Coulomb
interaction and raise `Tc` instead of lowering it.
"""
function checkMucsFromMu(inp, ef, omega_ph, phonon_cutoff)
    bad = Pair{Symbol,Float64}[]
    inp.muc_AD < 0 && push!(bad, :muc_AD => inp.muc_AD)
    inp.muc_ME < 0 && push!(bad, :muc_ME => inp.muc_ME)
    isempty(bad) && return nothing

    sig(x) = string(round(x, sigdigits=5))

    text = "The conversion from μ to μ* gave a negative pseudopotential:\n"
    for (name, val) in bad
        text *= "\n  * " * string(name) * " = " * sig(val)
    end
    text *= "\n\nμ* = μ / (1 + μ·ln(ε_el/ω)) is only meaningful while the denominator stays positive,"
    text *= "\ni.e. while the typical electron energy ε_el is well above the cutoff frequency ω."
    text *= "\nHere ε_el = " * sig(ef) * " meV against ω = " * sig(omega_ph) * " meV (μ*_AD)"
    text *= " and ω = " * sig(phonon_cutoff) * " meV (μ*_ME).\n"
    text *= "\nCheck, in this order:"
    text *= "\n  * typEl = " * sig(inp.typEl) * " meV - the typical electron energy. It must be in meV;"
    text *= "\n    a value left in eV or Ry is the most common cause. When typEl is unset, ef or efW is used."
    text *= "\n  * omega_c = " * sig(inp.omega_c) * " meV - enters μ*_ME; too large a cutoff drives the denominator negative."
    text *= "\n  * mu = " * sig(inp.mu) * " - the Coulomb strength N(ε_F)·W(ε_F,ε_F)."
    text *= "\n\nSetting muc_AD (and muc_ME) directly skips the conversion."
    text *= "\nSee the μ* conversion section of the Troubleshooting page and the pseudopotential section of the Input documentation.\n\n"

    error(text)
end
