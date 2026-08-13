"""
Format and write the output into 
    - log-file
    - Summary-file
    - Self Energy-file

Julia Packages:
    - 

Comments:
    -

"""

"""
    TableSpec(header, width, precision)

One console output table: column `header`s, column `width`s and `precision`s, plus the
lazily-filled horizontal rule `Hline` (set in `printTableHeader`) and `partingLine`
(set in `printStartMessage`).
"""
mutable struct TableSpec
    header::Vector{String}
    width::Vector{Int64}
    precision::Vector{Int64}
    Hline::String
    partingLine::String
end
TableSpec() = TableSpec(String[], Int64[], Int64[], "", "")
TableSpec(header, width, precision) = TableSpec(header, width, precision, "", "")

"""
    Console()

Console output state: a `cDOS` and a `vDOS` `TableSpec` (accessed uniformly for both
solvers; the imaginary-axis solver populates only the one selected by `cDOS_flag`) plus
the top-level `partingLine`.
"""
mutable struct Console
    cDOS::TableSpec
    vDOS::TableSpec
    partingLine::String
end
Console() = Console(TableSpec(), TableSpec(), "")

"""
    savePlotData(filepath, header, columns...)

Save tabular plot data as a tab-separated file. `header` is written as the first
line (prefix a `#` yourself if desired). Each element of `columns` is one data column.
"""
function savePlotData(filepath, header, columns...)
    open(filepath, "w") do io
        write(io, header * "\n")
        writedlm(io, zip(columns...), '\t')
    end
end


"""
    printAsciiArt()

print IsoME as Ascii art
"""
function printIsoME()
    strIsoME = "\n\n"
    strIsoME = strIsoME*"   _                 __  __   ______ \n"
    strIsoME = strIsoME*"  | |               |  \\/  | |  ____|\n"
    strIsoME = strIsoME*"  | |  ___    ___   | \\  / | | |__   \n"
    strIsoME = strIsoME*"  | | / __|  / _ \\  | |\\/| | |  __|  \n"
    strIsoME = strIsoME*"  | | \\__ \\ | (_) | | |  | | | |____ \n"
    strIsoME = strIsoME*"  |_| |___/  \\___/  |_|  |_| |______|\n\n"

    print(strIsoME)

    return strIsoME
end


"""
    printStartMessage(console)

Start message - Eliashberg Solver
"""
function printStartMessage(console::Console, inp, log_file; mode = 0)

    strAuthors =  "  Authors: Christoph Heil, Dominik Spath, Eva Kogler\n\n"

    # parting line of each sub-table (empty table -> trivial line, never used)
    console.cDOS.partingLine = "-"^(sum(console.cDOS.width) + length(console.cDOS.width) + 1)
    console.vDOS.partingLine = "-"^(sum(console.vDOS.width) + length(console.vDOS.width) + 1)

    # active table (cDOS_flag selects it for both solvers)
    strLine = inp.cDOS_flag == 1 ? console.cDOS.partingLine : console.vDOS.partingLine

    strMode = mode == 0 ? "Eliashberg Solver started" : "Real Axis Solver started"

    printTee(log_file, strAuthors)
    printTee(log_file, strLine)

    printTextCentered(strMode, strLine, file = log_file, bold = true)
    printTee(log_file, strLine*"\n\n\n")

    console.partingLine = strLine

    return console

end


"""
    printADtable(console)

Print the Allen-Dynes results to the console.
"""
function printADtable(console, ML_Tc, AD_Tc, BCS_gap, lambda, omega_log, log_file)

    # hline in rest of console
    partingLineCons = console.partingLine

    # table headline
    headline = "Allen-Dynes-McMillan Formula";
    delimiter = "|"

    # table Results
    ADvalues = [(round(max(ML_Tc, 0), digits=2)) (round(max(AD_Tc, 0), digits=2)) (round(BCS_gap, digits=2)) (round(lambda, digits=2)) (round(omega_log, digits=2))]

    # table header
    header  = ["Tc_ML", "Tc_AD", "BCS_gap", "lambda", "omega_log"]
    units = [" K", " K", " meV", "", " meV"]

    # calc width of each column
    width = Vector{Int}(zeros(length(header)))
    for k in eachindex(header)
        width[k] = maximum([length(header[k])+2, length(string(ADvalues[k])) + length(units[k]) + 2])
    end
    header = formatTableHeader(header, width)

    # Hline
    Hline = "."*"-"^(length(width)-1)
    for w in width
        Hline = Hline*"-"^w
    end
    Hline= Hline*"."
    blanksAD = Int(maximum([0, floor((length(partingLineCons)-length(Hline) )/2)]))

    # format headline
    lenLeft = Int(ceil((length(Hline) - length(headline)-2)/2))
    lenRight = Int(floor((length(Hline) - length(headline)-2)/2))
    headline = " "^lenLeft*headline*" "^lenRight

    # print headline
    logText = " "^blanksAD*Hline*"\n"*" "^blanksAD*delimiter*headline*delimiter*"\n"
    println(" "^blanksAD*Hline)
    print(" "^blanksAD*delimiter)
    printstyled(headline; bold=true)
    println(delimiter)
   
    # Hline
    logText = logText*" "^blanksAD*replace(Hline, "." => "|")*"\n"
    println(" "^blanksAD*replace(Hline, "." => "|"))

    # Define table header 
    logText = logText*" "^blanksAD*delimiter
    print(" "^blanksAD*delimiter)
    for k in eachindex(header)
        value = header[k]
        # print
        printstyled(value; bold=true)
        print(delimiter)

        # save for log file
        logText = logText*string(value)*delimiter
    end

    ### parting line ###
    logText = logText*"\n"*" "^blanksAD*replace(Hline, "." => "|")*"\n"*" "^blanksAD
    print("\n"*" "^blanksAD*replace(Hline, "." => "|")*"\n"*" "^blanksAD)

    ### AD values ###
    for k in eachindex(ADvalues)
        value = ADvalues[k]
        w = width[k]
        unit = units[k]

        numDig = numDigits(value)
        blanks = (w - (sum(numDig) + 1) - length(unit)) / 2

        print(delimiter)
        print(" "^Int(floor(blanks)), string(value), unit, " "^Int(ceil(blanks)))
        logText = logText*"|"*" "^Int(floor(blanks))* string(value)* unit* " "^Int(ceil(blanks))
    end
    logText = logText*"|\n"*" "^blanksAD*replace(Hline, "." => " ")*"\n\n"
    print("|\n"*" "^blanksAD*replace(Hline, "." => " ")*"\n\n")

    # write everything to log file
    print(log_file, logText*"\n")
end


"""
    printSummary()

Summarize Tc calculation and print it to the console
"""
function printSummary(inp, Tc, log_file)

    text = ""
    if Tc[2] <= 0.5
        text = text * "\n - " * inp.material * " is not a superconductor above T = 0.5 K"
    elseif isnan(Tc[1])
        text = text * "\n - Couldn't find a superconducting gap in the specified area"
        text = text * "\n - Consider searching below T = " * string(Tc[2]) * " K"
    elseif isnan(Tc[2])
        text = text * "\n - " * inp.material * " is a superconductor"
        text = text * "\n - Highest given temperature reached"
        text = text * "\n - Tc > " * string(Tc[1]) * " K"
    else
        text = text * "\n - " * inp.material * " is a superconductor"
        text = text * "\n - Tc = " * string(round((Tc[2]+Tc[1])/2, digits=2)) * " (±"* string(round((Tc[2]-Tc[1])/2, digits=2)) *")" * " K"
    end

    printstyled("\nSummary:", bold=true)
    println(text*"\n")

    # log file
    print(log_file, "\nSummary:"*text*"\n")

end


"""
    printTextCentered(text, hline[, boldFlag, blanks, delimiter])

print a text centered within a line consisting of delimiters
"""
function printTextCentered(text, hline; file = "", bold = false, blanks=3, delimiter = "-", newline = "\n", consoleFlag = true)

    lenLeft = Int(ceil((length(hline) - length(text))/2) - blanks)
    lenRight = Int(floor((length(hline) - length(text))/2) - blanks)

    leftText = newline*delimiter^lenLeft*" "^blanks
    rightText = " "^blanks * delimiter^lenRight*"\n"

    # print to console
    if consoleFlag
        print(leftText)
        if bold
            printstyled(text; bold=true)
        else
            print(text)
        end
        print(rightText)
    end

    # print to file
    if ~isa(file, IOBuffer) && isfile(file)
        print(file, leftText)
        print(file, text)
        print(file, rightText)
    end

end


"""
    printTableHeader(table, initValues, log_file)

Initialize the table header (setting `table.Hline`) and print it, together with the
initial-value row `initValues`, to the console and log file.
"""
function printTableHeader(table::TableSpec, initValues, log_file)

    width = table.width
    header = table.header

    ### Define boundary ###
    tableHline = ""
    for w in width
        tableHline = tableHline*"."*"-"^w
    end
    tableHline = tableHline*"."
    println(tableHline)

    ### Define table header ###
    delimiter = "|"
    print(delimiter)
    tableHeader = tableHline*"\n"*delimiter
    for k in eachindex(header)
        value = header[k]
        # print
        printstyled(value; bold=true)
        print(delimiter)

        # save for log file
        tableHeader = tableHeader*string(value)*delimiter
    end

    ### parting line ###
    tableHeader = tableHeader*"\n"*replace(tableHline, "." => "|")
    print("\n"*replace(tableHline, "." => "|")*"\n")

    ### Initial values ###
    initValues, strFormat, format = formatTableRow(initValues, width, 2)
    printTableRow(stdout, initValues, strFormat, format)

    # write everything to log_file
    print(log_file, tableHeader * "\n")
    printTableRow(log_file, initValues, strFormat, format)

    table.Hline = tableHline

    return table

end


"""
    formatTableHeader(console)

Format the header of the console output s.t. each column has the 
specified length
"""
function formatTableHeader(console::Console)
    formatTableHeader(console.cDOS)
    formatTableHeader(console.vDOS)
    return console
end

function formatTableHeader(t::TableSpec)
    for k in eachindex(t.header)
        blanks = (t.width[k] - length(t.header[k])) / 2
        t.header[k] = " "^Int(floor(blanks)) * t.header[k] * " "^Int(ceil(blanks))
    end
    return t
end

function formatTableHeader(header, width)

    for k in eachindex(header)
        blanks = (width[k]-length(header[k]))/2
        header[k] = " "^Int(floor(blanks))*header[k]*" "^Int(ceil(blanks))
    end

    return header
end


"""
    formatTableRow(vec, widthCol, prec=5, logConsole=true)
    
Format the console output s.t. it is aligned to the header
"""
function formatTableRow(vec, widthCol, prec=5, logConsole=true)
    """
    Format the console output s.t. it is aligned to the header
    This is done by introducing a dynamic width and precision via
    "%*s" (requires Julia 1.10)
    The output can be printed via printf   
    Remark: The use of printf could be circumvented by using just " "^width

    -------------------------------------------------------------------
    Input:
        vec:        vector, each element should be written to the console
        widthCol:   available length of each column
        prec:       precision after the comma, optional
        strAttach:  attach this string to vec

    --------------------------------------------------------------------
    Output:
        vec:        rounded to precision    
        out:        formatting of vec given as string suitable for 
                    printf function, e.g "%*s"
        format:     whitespace formating for out

    --------------------------------------------------------------------
    Comments:
        - Use of dynamic width and precision requires Julia 1.10
        - To print the output using a dynamic width and precision use
            for i in axes(strConsole, 1)
                Printf.format(stdout, Printf.Format(out[i]), format[i, 1], " ", format[i, 2], format[i, 3], vec[i], format[i, 4], " ")
            end

    --------------------------------------------------------------------
    """


    if length(prec) == 1
        prec = fill(prec[1], length(vec))
    end

    out = Array{String}(undef, length(vec))
    format = zeros(Int, length(vec), 4)
    for k in eachindex(vec)
        value = vec[k]
        width = widthCol[k]

        if isnothing(value) || !isfinite(value)
            # NaN/Inf get a "-" placeholder; numDigits() cannot parse them either
            numDig = [0,0]
            spacing = (width - (sum(numDig) + 1)) / 2
            format[k, :] = [Int(floor(spacing)), Int(numDig[1]), Int(numDig[2]), Int(ceil(spacing))]
            out[k] = " "^Int(format[k,1])*"-"*" "^Int(format[k,4])*"|"
        
        else
            value = round(value, digits=prec[k])
            numDig = numDigits(value)
            if (numDig[1] + prec[k]) > width
                numDig[2] = minimum([0, width - numDig[1]])
            else
                numDig[2] = prec[k]
            end

            if k == 1    # nr iteration, no comma
                out[k] = "|%*s%*.*f%*s|"
                if logConsole
                    spacing = (width - (numDig[1])) / 2
                    format[k, :] = [Int(floor(spacing)), Int(numDig[1]), Int(0), Int(ceil(spacing))]
                else
                    spacing = (width - (sum(numDig) + 1)) / 2
                    format[k, :] = [Int(floor(spacing)), Int(numDig[1]), Int(numDig[2]), Int(ceil(spacing))]
                end
            else
                out[k] = "%*s%*.*f%*s|"
                spacing = (width - (sum(numDig) + 1)) / 2
                format[k, :] = [Int(floor(spacing)), Int(numDig[1]), Int(numDig[2]), Int(ceil(spacing))]
            end
        end
    end

    out[end] = out[end]*"\n"

    return vec, out, format        

end


"""
    printTableRow(io, vec, strFormat, format)

Write one formatted table row (the output of [`formatTableRow`](@ref)) to `io`.

Entries that are `nothing` or non-finite are rendered by `formatTableRow` as a plain "-"
placeholder, i.e. a string *without* any format specifier. Those must be printed verbatim -
handing them to `Printf.format` together with the seven width/precision arguments throws
`ArgumentError: Number of format specifiers and number of provided args differ: 0 != 7`.
"""
function printTableRow(io::IO, vec, strFormat, format)
    for i in axes(strFormat, 1)
        value = vec[i]
        if isnothing(value) || !isfinite(value)
            print(io, strFormat[i])
        else
            Printf.format(io, Printf.Format(strFormat[i]), format[i, 1], " ", format[i, 2], format[i, 3], value, format[i, 4], " ")
        end
    end
    return nothing
end


"""
    numDigits(x)

Gives the number of digits before and after the comma
"""
function numDigits(x)
    digits = zeros(Int, length(x), 2)
    for k in eachindex(x)
        str = split(string(x[k]), ".")
        digits[k,:] = [length(str[1]), length(str[2])]
    end
    return digits

end


"""
    printFlagsAsText()

Print the flag values as text to the console
"""
function printFlagsAsText(inp, log_file; mode ="Matsubara")
    text = ""
    if inp.material != "Material"
        text *=  " - Material: "*inp.material*" \n"
    end

    # search mode
    if inp.temps == [-1]
        text *= " - Tc search mode activated\n"
    else 
        if length(inp.temps) == 1
            text *= " - Tc search range: "*string(inp.temps[1])*" K\n"
        else
            text *= " - Tc search range: "*string(minimum(inp.temps))*" - "*string(maximum(inp.temps))*" K\n"
        end

    end
    
    # cut off
    if mode == "Matsubara"
        text *= " - Matsubara cutoff: "*string(inp.imOmega_c)*" meV\n"
    elseif mode == "realFreq"
        text *= " - Frequency cutoff: "*string(inp.reOmega_c)*" meV\n"
        text *= " - Frequency grid step: "*string(inp.domega)*" meV\n"
    end

    # cDos
    if inp.cDOS_flag == 0
        if inp.mu_flag == 1
            text *= " - Variable DoS with μ-update\n"
        else
            text *= " - Variable DoS with constant μ = ϵ_F\n"
        end
    elseif inp.cDOS_flag == 1
        text *= " - Constant DoS\n"
    end

    # Weep
    if inp.include_Weep == 1
        text *= " - Static Coulomb interaction W(e,ep) in "*inp.Weep_unit*"\n"
    elseif inp.include_Weep == 0
        text *= " - Morel-Anderson pseudopotential\n"
        text *= "     - μ*_AD = "*string(round(inp.muc_AD, digits=3))*"\n"
        text *= "     - μ*_ME = "*string(round(inp.muc_ME, digits=3))*"\n"
    end

    text *= "\n"

    print(text)
    
    # log file
    print(log_file, text)


end


"""
"""
function printTee(log_file, text)
    print(text)
    print(log_file, text)
    flushLog(log_file)
end


"""
    flushLog(log_file)

Flush the log stream. The log file is buffered, so a run that is killed from the
outside (OOM killer, SIGTERM, segfault in a library) would otherwise lose
everything written so far. Never throws.
"""
function flushLog(log_file)
    try
        isopen(log_file) && flush(log_file)
    catch
        # a failing flush must never mask the reason we are flushing
    end
    return nothing
end


"""
    closeLog(log_file)

Flush and close the log stream if it is still open. Never throws.
"""
function closeLog(log_file)
    flushLog(log_file)
    try
        isopen(log_file) && close(log_file)
    catch
    end
    return nothing
end

"""
    printError(text, ex, log_file, errorLogger)

print a formatted error message to the console and log_file
"""
function printError(text, ex, log_file, errorLogger)

    print(log_file, "\n")
    with_logger(errorLogger) do
        @error text exception = ex
    end
    print(log_file, "\nFor further information please refer to the CRASH file\n\n")
    closeLog(log_file)

    print("\n")
    rethrow(ex)
end


"""
    printWarning(text, ex, log_file)

print a formatted warning message to the console and log_file
"""
function printWarning(text, log_file; ex = nothing)

    printTee(log_file, "\n")
    if isnothing(ex)
        @warn text
        printTee(log_file, "\n")
    else
        @warn text exception = ex
        printTee(log_file, "For further information please refer to the CRASH file\n")
    end

end


"""
    IsoMEError(what, cause, bt)

A failure in one of the solver stages, tagged with the message that tells the
user *where* it happened ("while reading the inputs", ...). It carries the
original exception and its backtrace up to the global error handler, which is
the single place that reports it. Marking the stage this way needs no shared
state, so it stays correct if the solver is ever parallelised.
"""
struct IsoMEError <: Exception
    what::String
    cause::Any
    bt::Any
end


"""
    stage(f, what)

Run a mandatory step of the solver. If it fails, the exception is tagged with
`what` and passed on to the global error handler, which writes the CRASH file and
the error message. Does no I/O itself.

    stage("in input structure") do
        checkInput!(inp)
    end
"""
function stage(f::F, what::AbstractString) where {F}
    try
        return f()
    catch ex
        ex isa InterruptException && rethrow(ex)     # Ctrl+C is not a solver failure
        ex isa IsoMEError && rethrow(ex)             # keep the innermost, most specific message
        throw(IsoMEError(what, ex, catch_backtrace()))
    end
end


"""
    attempt(f, inp, log_file, what)

Run an optional step of the solver. If it fails, the exception is written to the
CRASH file and reported as a warning, and the run continues. Reports on the spot
rather than deferring, because execution must not leave this point.

    attempt(inp, log_file, "Error while plotting. Skipping plots.") do
        createFigures(inp, matval, Delta0, temps, Tc, log_file)
    end
"""
function attempt(f::F, inp, log_file, what::AbstractString) where {F}
    try
        f()
    catch ex
        ex isa InterruptException && rethrow(ex)     # Ctrl+C must stop the run, not be skipped
        writeToCrashFile(inp, ex, catch_backtrace())
        printWarning(what * " This step is optional and the run continued; the results are still valid. If it keeps failing, please open an issue at https://github.com/cheil/IsoME.jl/issues and attach the CRASH file.", log_file, ex = ex)
    end
    return nothing
end


"""
    handleFatalError(ex, bt, inp, log_file, errorLogger)

Global error handler: the single place where a fatal error is written to the
CRASH file and to the log. Handles both errors tagged by `stage` and errors that
escaped the stages altogether, e.g. from the solver driver itself. Never throws -
the caller decides what to rethrow.
"""
function handleFatalError(ex, bt, inp, log_file, errorLogger)

    # a deliberate abort is not a crash: only a note, no error report
    if ex isa InterruptException
        writeNoteToCrashFile(inp, "User interruption")
        printTee(log_file, "\nUser interruption\n")
        return nothing
    end

    if ex isa IsoMEError
        text  = ex.what * ". Stopping now!"
        cause = ex.cause
        trace = ex.bt                   # where it actually broke, not the unwound view
    else
        text  = "Unexpected internal error. Stopping now!"
        cause = ex
        trace = bt
    end

    writeToCrashFile(inp, cause, trace)

    try
        print(log_file, "\n")
        with_logger(errorLogger) do
            # only the exception, not the backtrace - the full trace is in the
            # CRASH file, and it would flood the log file
            @error text exception = cause
        end
        print(log_file, "\nFor further information please refer to the CRASH file\n\n")
    catch
        # the log file itself is broken - the CRASH file written above is what is left
    end
    flushLog(log_file)

    print("\n")
    return nothing
end


"""
    crashFilePath(inp)

Path of the CRASH file. Falls back to the working directory if `outdir` does not
exist yet, so that a crash during the setup is not silently lost.
"""
function crashFilePath(inp)
    dir = isempty(inp.outdir) ? pwd() : inp.outdir
    isdir(dir) || (dir = pwd())
    return joinpath(dir, "CRASH")
end


"""
    writeToCrashFile(inp, ex, bt)

Save exception and backtrace in CRASH file. Never throws: if the CRASH file can
not be written the report goes to stderr instead, so that a failure here can not
mask the original error.
"""
function writeToCrashFile(inp, ex, bt)
    crashFile = nothing
    try
        crashFile = inp.testMode ? IOBuffer() : open(crashFilePath(inp), "a")
        println(crashFile, "="^80)
        showerror(crashFile, ex, bt)
        print(crashFile, "\n\n")
        flush(crashFile)
    catch
        try
            print(stderr, "\n[IsoME] Could not write the CRASH file. Original error:\n")
            showerror(stderr, ex, bt)
            print(stderr, "\n\n")
        catch
        end
    finally
        crashFile === nothing || close(crashFile)
    end
    return nothing
end


"""
    writeToCrashFile(inp)

Save the exception that is currently being handled in CRASH file. Only valid
inside a catch block. Never throws.
"""
function writeToCrashFile(inp)
    crashFile = nothing
    try
        crashFile = inp.testMode ? IOBuffer() : open(crashFilePath(inp), "a")
        print(crashFile, current_exceptions())
        print(crashFile, "\n\n")
        flush(crashFile)
    catch
        print(stderr, "\n[IsoME] Could not write the CRASH file.\n")
    finally
        crashFile === nothing || close(crashFile)
    end
    return nothing
end


"""
    writeNoteToCrashFile(inp, text)

Save a plain note in the CRASH file, for the cases that end the run without being
an exception to report - a Ctrl+C in particular. Never throws.
"""
function writeNoteToCrashFile(inp, text)
    crashFile = nothing
    try
        crashFile = inp.testMode ? IOBuffer() : open(crashFilePath(inp), "a")
        println(crashFile, "="^80)
        println(crashFile, text)
        print(crashFile, "\n\n")
        flush(crashFile)
    catch
        print(stderr, "\n[IsoME] Could not write the CRASH file.\n")
    finally
        crashFile === nothing || close(crashFile)
    end
    return nothing
end


"""
    writeInputFlags(Tc ,inp, out_vars, header)

Save results of each iteration and input parameters in a file Info.txt
"""
function createInfoFile(inp)
    # write to output file
    name = "Info.txt"

    if isfile(inp.outdir*name)
        rm(inp.outdir*name)
    end
    outfile = open(inp.outdir*name, "w")
    
    ### Input parameters ###
    print(outfile, replace(replace(join(inp.all, "\n"), "nothing"=>"-"), "[-1]"=>"-"))

    close(outfile)
end


"""
    summarizeResults(Tc, out_vars, header)

Save Delta(0) at each temperature
"""
function createSummaryFile(inp::arguments, Tc, out_vars, header)
    # write to summary file
    name = "Summary.dat"

    if isfile(inp.outdir*name)
        rm(inp.outdir*name)
    end
    outfile = open(inp.outdir*name, "w")

    ### header
    out = "# "
    if isnan(Tc[1])
        out = out * "Tc < " * string(Tc[2]) *" K"
    elseif isnan(Tc[2])
        out = out * "Tc > " * string(Tc[1]) * " K"
    else
        out = out * "Tc = " * string(round((Tc[2]+Tc[1])/2, digits=2)) * " (±"* string(round((Tc[2]-Tc[1])/2, digits=2)) *")" * " K"
    end
    out = out*"\n"*header*"\n"
    print(outfile, out)

    close(outfile)
      
    # save header & gap
    #writedlm(outfile, round.(out_vars, digits=4), '\t')
    spacingsHeader = findall("   ", header)
    open(inp.outdir*name, "a") do io
        for row in eachrow(out_vars)
            startCol = 0
            for (val, spacing) in zip(row, spacingsHeader)
                space = spacing[end] - startCol
                @printf(io, "%-*.3f", space, val)       # could be changed to print whole row at once + less digits after comma for T 
                startCol = spacing[end]+1
            end
            println(io)
        end
    end

end



function createFigures(inp, matval, Delta0, temps, Tc, log_file)

    # values
    a2f_omega_fine, a2f_fine = matval

    # defaults
    plot_font = "Computer Modern"
    default(
        fontfamily=plot_font,
        linewidth=2,
        framestyle=:box,
        label=nothing,
        grid=false
    )

    # print a2F vs. energy
    xlim_max = Int(round(maximum(a2f_omega_fine) / 10 * 1.01, RoundUp) * 10)
    xtick_val = 0:10:xlim_max
    ylim_max = Int(round(maximum(a2f_fine), RoundUp))

    plot(a2f_omega_fine, a2f_fine,1)
    xlims!(0, xlim_max)
    ylims!(0, ylim_max)
    if inp.material != "Material"
        title!(inp.material)
    end
    xlabel!(L"\omega ~ \mathrm{(meV)}")
    ylabel!(L"\alpha^2F ~ \mathrm{(1)}")
    savefig(inp.outdir * "/a2F_sm" * string(inp.ind_smear) * ".pdf")

    if all(isnan.(Delta0))
        printstyled("Info: "; color=:blue)
        println("No superconducting gap found - skipping plot\n")

        print(log_file,  "Info: No superconducting gap found - skipping plot\n")
    else
        # print gap vs. temperature. Each of these is assigned exactly once: `temps_plot` is
        # captured by the fit model below, and reassigning a captured variable boxes it,
        # which would leave the whole fit and every axis limit untyped
        keep = .~isnan.(Delta0)
        order = sortperm(temps[keep])
        temps_plot = temps[keep][order]
        Delta0_plot = Delta0[keep][order]

        if maximum(temps_plot) < 10
            xlim_max = round(maximum(temps_plot) * 1.1, RoundUp)
            xtick_val = 0:1:xlim_max
        elseif maximum(temps_plot) < 20
            xlim_max = round(maximum(temps_plot) * 1.1, RoundUp)
            xtick_val = 0:2:xlim_max
        else
            xlim_max = round(maximum(temps_plot) / 10 * 1.01, RoundUp) * 10
            xtick_val = 0:10:xlim_max
        end
        ylim_max = round(maximum(Delta0_plot)*1.11, RoundUp)
        
        # Create a custom gradient
        my_gradient = cgrad(:coolwarm)  
        #gradientBlue = cgrad(:Blues, rev=true)
        #gradientRed = cgrad(:Reds)
        # Define gradient range
        max_gradient_val = 77               # Blue Gradient applies up to this value
        
        # normalize x: blue at 0 K, white at max_gradient_val, red at 2*max_gradient_val.
        # If the gap survives beyond 2*max_gradient_val the range is stretched to the
        # highest temperature instead, so the colour scale never wraps back to blue.
        color_values = temps_plot ./ max(2 * max_gradient_val, maximum(temps_plot))

        # Map colors: Use blue gradient for values <= max_gradient_val, red gradient for rest
        marker_colors = [my_gradient[clamp(v, 0.0, 1.0)] for v in color_values]
        
        h=scatter(temps_plot, Delta0_plot, color=marker_colors, colorbar=false, markerstrokewidth=1, ms=6, xticks=xtick_val)
        if length(temps_plot) > 1
            try
                #p[1] = Delta[1], exp(p[2])+max(temps_plot) = Tc, p[3] = fit parameter adjusting the curvature
                Delta(T,p) = p[1]* tanh.((π*kb*(exp(p[2])+maximum(temps_plot)))/p[1]*sqrt.(p[3] * ((exp(p[2])+maximum(temps_plot))./T .− 1) ))
                pGuess = convert(Vector{Float64}, [(Delta0_plot[1]), 0, 1])
            
                fit = curve_fit(Delta, temps_plot, Delta0_plot, pGuess)
                par = fit.param
                TcFit = exp(par[2])+maximum(temps_plot)
                Delta0Fit = abs(par[1])     # Delta(T) is symmetric wrt Delta0

                if ~fit.converged
                    printWarning("The Δ(T) fit did not converge within $(fit.iterations) iterations "*
                                 "(Tc_fit = $(round(TcFit, digits=2)) K, Δ₀_fit = $(round(Delta0Fit, digits=3)) meV). "*
                                 "The fitted curve may be unreliable.", log_file)
                end

                # second derivative
                A = π*kb*TcFit/Delta0Fit
                a = par[3]
                DeltaSecDer(T) = @. -(a^2 * A * TcFit * sech(A * sqrt(a * (TcFit/T - 1)))^2 * 
                        (-3*TcFit + 4*T + 2*A * TcFit * sqrt(a * (TcFit - T) / T) * 
                        tanh(A * sqrt(a * (TcFit/T - 1))))) / 
                        (4 * (a * (TcFit - T) / T)^(3/2) * T^4)

                # 
                temps_fit = collect(range(0.01, TcFit, 100))
                Delta_fit = Delta(temps_fit, par)
                
                # plot if negative curvature
                negCurv = all(DeltaSecDer(temps_fit).<0)
                aboveTc = TcFit > Tc[1]
                if negCurv && aboveTc #&& (TcFit < Tc[2] || isnan(Tc[2])) && TcFit < 1e3 # plot only if TcFit is reasonable ??
                    plot!(temps_fit, Delta_fit, linestyle=:dash, linewidth=1, color=:gray, z_order=:back)
            
                    # adjust limits
                    xlim_max = maximum([xlim_max, TcFit+1])
                    ylim_max = maximum([ylim_max, Delta0Fit+1])
            
                    if xlim_max <= 10
                        xtick_val = 0:1:xlim_max
                    elseif xlim_max > 1e3
                        error("TcFit unreasonable")
                    elseif xlim_max <= 20
                        xtick_val = 0:2:xlim_max
                    else
                        xtick_val = 0:10:xlim_max
                    end
                else
                    # fit returned, but the result is not physical enough to show
                    reason = String[]
                    negCurv || push!(reason, "the curvature of Δ(T) is not negative over the whole fit range")
                    aboveTc || push!(reason, "the fitted Tc ($(round(TcFit, digits=2)) K) does not exceed the "*
                                             "lower bound of the Tc interval ($(round(Tc[1], digits=2)) K)")
                    printWarning("Not plotting the Δ(T) fit because "*join(reason, " and ")*
                                 ". The calculated gap values are shown without a fit curve.", log_file)
                end
            catch  ex
                ex isa InterruptException && rethrow(ex)
                # handled here rather than in the enclosing attempt(), so the rest of the
                # figure is still produced - but keep attempt()'s CRASH file, otherwise the
                # warning would point at a file that was never written
                writeToCrashFile(inp, ex, catch_backtrace())
                printWarning("The Δ(T) fit failed. The calculated gap values are shown without a fit curve.",
                             log_file, ex=ex)
            end
        end
        plot!(xticks=xtick_val)
        xlims!(0, xlim_max)
        ylims!(0, ylim_max)
        xlabel!(L"T ~ [\mathrm{K}]")
        ylabel!(L"\Delta_0 ~ [\mathrm{meV}]")
        
        #display(h)
        
        namePlot = "Delta0"
        if inp.material != "Material"
            namePlot = namePlot * "_" * inp.material
        end
        if inp.cDOS_flag == 1
            namePlot = namePlot * "_" * "cDOS"
        else
            namePlot = namePlot * "_" * "vDOS"
        end
        if inp.include_Weep == 1
            namePlot = namePlot * "_" * "W"
        else
            namePlot = namePlot * "_" * "muc"
        end
        namePlot = namePlot * ".pdf"
        savefig(inp.outdir * namePlot)
    end
end



"""
    writeSelfEnergy(itemp, inp, components)

Write each self-energy component to `outdir/SelfEnergy/` as both a `.dat` file and a `.png` plot.
Single entry point for both axes and all modes; every component carries its own grid, column header
and axis labels, so the frequency label (ω vs iωₙ) and the units are always correct.

`components` is a vector of NamedTuples `(name, x, y, xlabel, ylabel, xhead, yhead)`. A complex `y` is
written as Re/Im columns and plotted as Re (blue) / Im (red); a real `y` (imaginary-axis quantities)
gets a single column and curve.
"""
function writeSelfEnergy(itemp, inp, components)

    folder = inp.outdir * "SelfEnergy/"
    if ~isdir(folder)
        mkdir(folder)
    end

    T = Float64(itemp)                          # consistent file names

    # `components` is heterogeneous, so its element type is the abstract NamedTuple. The
    # per-component work goes through a function barrier: each element is dispatched once
    # and writeSelfEnergyComponent then runs on its own concrete NamedTuple type
    for c in components
        writeSelfEnergyComponent(folder, T, inp.material, c)
    end

    return nothing
end


"""
    writeSelfEnergyComponent(folder, T, material, c)

Write one self-energy component of [`writeSelfEnergy`](@ref) as `.dat` + `.png`.
"""
function writeSelfEnergyComponent(folder::AbstractString, T::Float64, material::AbstractString, c::NamedTuple)

    base = folder * c.name * "_" * string(T) * "K"

    # ----- data file -----
    open(base * ".dat", "w") do io
        if eltype(c.y) <: Complex
            write(io, "#  " * c.xhead * "\tRe(" * c.yhead * ")\tIm(" * c.yhead * ")\n")
            writedlm(io, [c.x real.(c.y) imag.(c.y)], '\t')
        else
            write(io, "#  " * c.xhead * "\t" * c.yhead * "\n")
            writedlm(io, [c.x c.y], '\t')
        end
    end

    # ----- plot -----
    if eltype(c.y) <: Complex
        plot(c.x, real.(c.y), color = :blue, label = "Real", linewidth = 2, xlabel = c.xlabel, ylabel = c.ylabel)
        plot!(c.x, imag.(c.y), color = :red, label = "Imag", linewidth = 2)
    else
        plot(c.x, c.y, color = :blue, label = "", linewidth = 2, xlabel = c.xlabel, ylabel = c.ylabel)
    end
    material != "Material" && title!(material)
    savefig(base * ".png")

    return nothing
end


"""
    saveSelfEnergy_matsubara(itemp, inp, wsi, deltai, znormi; chi, phiph, phic, epsilon)

Assemble the imaginary-axis (Matsubara) self-energy components and hand them to
[`writeSelfEnergy`](@ref). Covers cDOS/vDOS and the μ / W(ε,ε′) modes: χ, φ_ph, φ_c are written only
when provided. In vDOS+W, `deltai` is Δ(ε_F, iωₙ) — the full ε-dependence is recoverable from the
φ_ph(iωₙ) and φ_c(ε) files.
"""
function saveSelfEnergy_matsubara(itemp, inp, wsi, deltai, znormi; chi = nothing, phiph = nothing, phic = nothing, epsilon = nothing)

    comps = NamedTuple[]

    push!(comps, (name = "Z", x = wsi, y = znormi,
                  xlabel = "iωₙ / meV", ylabel = "Z(iωₙ) / 1", xhead = "iωₙ / meV", yhead = "Z(iωₙ) / 1"))

    # Δ carries an ε-dependence exactly where φ_c does; the caller passes the ε_F slice
    deltaLabel = isnothing(phic) ? "Δ(iωₙ) / meV" : "Δ(ε_F, iωₙ) / meV"
    push!(comps, (name = "Delta", x = wsi, y = deltai,
                  xlabel = "iωₙ / meV", ylabel = deltaLabel, xhead = "iωₙ / meV", yhead = deltaLabel))

    isnothing(chi)   || push!(comps, (name = "Chi", x = wsi, y = chi,
                  xlabel = "iωₙ / meV", ylabel = "χ(iωₙ) / meV", xhead = "iωₙ / meV", yhead = "χ(iωₙ) / meV"))
    isnothing(phiph) || push!(comps, (name = "Phiph", x = wsi, y = phiph,
                  xlabel = "iωₙ / meV", ylabel = "φ_ph(iωₙ) / meV", xhead = "iωₙ / meV", yhead = "φ_ph(iωₙ) / meV"))
    isnothing(phic)  || push!(comps, (name = "Phic", x = epsilon, y = phic,
                  xlabel = "ε / meV", ylabel = "φ_c(ε) / meV", xhead = "ε / meV", yhead = "φ_c(ε) / meV"))

    writeSelfEnergy(itemp, inp, comps)
    return nothing
end


"""
    saveSelfEnergy_realAxis(itemp, inp, w_static, delta, Z; chi, phi_ph, phi_c, epsilon)

Assemble the real-axis self-energy components (Δ, Z, χ and φ_ph all on `w_static`; φ_c(ε) on the
ε-grid) and hand them to [`writeSelfEnergy`](@ref). Optional channels are written only when
provided, so the same call serves cDOS and vDOS (+μ / +W).
"""
function saveSelfEnergy_realAxis(itemp, inp, w_static, delta, Z; chi = nothing, phi_ph = nothing, phi_c = nothing, epsilon = nothing)

    comps = NamedTuple[]

    push!(comps, (name = "Delta", x = w_static, y = delta,
                  xlabel = "ω / meV", ylabel = "Δ(ω) / meV", xhead = "ω / meV", yhead = "Δ(ω) / meV"))
    push!(comps, (name = "Z", x = w_static, y = Z,
                  xlabel = "ω / meV", ylabel = "Z(ω) / 1", xhead = "ω / meV", yhead = "Z(ω) / 1"))

    isnothing(chi) || push!(comps, (name = "Chi", x = w_static, y = chi,
                  xlabel = "ω / meV", ylabel = "χ(ω) / meV", xhead = "ω / meV", yhead = "χ(ω) / meV"))
    isnothing(phi_ph) || push!(comps, (name = "Phiph", x = w_static, y = phi_ph,
                  xlabel = "ω / meV", ylabel = "φ_ph(ω) / meV", xhead = "ω / meV", yhead = "φ_ph(ω) / meV"))
    (isnothing(phi_c) || isnothing(epsilon)) || push!(comps, (name = "Phic", x = epsilon, y = phi_c,
                  xlabel = "ε / meV", ylabel = "φ_c(ε) / meV", xhead = "ε / meV", yhead = "φ_c(ε) / meV"))

    writeSelfEnergy(itemp, inp, comps)
    return nothing
end



"""
     Base.getproperty(a::arguments, v::Symbol)

return all arguments in input structure arguments
"""
function Base.getproperty(a::arguments, v::Symbol)
    if v == :all
        input = Vector{String}()
        for name in fieldnames(arguments)
            text = string(name)*": "*string(getfield(a, name))
            input = push!(input, text)
        end

        return input
    elseif v == :args
        input = Vector{String}()
        for name in fieldnames(arguments)
            input = push!(input, string(name))
        end

        return input
    else
        return getfield(a, v)
    end
end


"""
     Base.setproperty!(a::arguments, v::Symbol, x)

set a field of the input structure, reporting unknown field names and wrong data types
in the same way as the `arguments` constructor
"""
function Base.setproperty!(a::arguments, v::Symbol, @nospecialize(x))
    hasfield(arguments, v) || _throw_unknown(arguments, [v])
    Ft = fieldtype(arguments, v)
    x isa Ft && return setfield!(a, v, x)
    try
        return setfield!(a, v, convert(Ft, x))
    catch ex
        ex isa InterruptException && rethrow(ex)
        _throw_invalid(arguments, Symbol[], Pair{Symbol,Any}[v => x])
    end
end
