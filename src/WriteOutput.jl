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
    printStartMessage(console)

Start message
"""
function printStartMessage(inp, console)

    strIsoMe = "\n\n"
    strIsoMe = strIsoMe*"   _                 __  __   ______ \n"
    strIsoMe = strIsoMe*"  (_)               |  \\/  | |  ____|\n"
    strIsoMe = strIsoMe*"   _   ___    ___   | \\  / | | |__   \n"
    strIsoMe = strIsoMe*"  | | / __|  / _ \\  | |\\/| | |  __|  \n"
    strIsoMe = strIsoMe*"  | | \\__ \\ | (_) | | |  | | | |____ \n"
    strIsoMe = strIsoMe*"  |_| |___/  \\___/  |_|  |_| |______|\n\n"

    strLine = "-"^(sum(console["width"])+length(console["width"])+1)

    strAuthors = "  Authors: Christoph Heil, Eva Kogler, Dominik Spath\n\n"

    strEliash = "Eliashberg Solver started"

    print(strIsoMe)
    #printTextCentered(inp, "Version 1.0", strLine, false)
    print(strAuthors)
    print(strLine)
    if inp.flag_log == 1
        print(inp.log_file, strIsoMe)
        print(inp.log_file, strAuthors)
        print(inp.log_file, strLine)
    end
    printTextCentered(inp, strEliash, strLine, true)
    print(strLine*"\n\n\n")
    if inp.flag_log == 1
        print(inp.log_file, strLine*"\n\n\n")
    end

    console["partingLine"] = strLine
    
    return console
                                    
end


"""
    printADtable(console)

Print the Allen-Dynes results to the console.
"""
function printADtable(inp, console, ML_Tc, AD_Tc, BCS_gap, lambda, omega_log)

    # hline in rest of console
    partingLineCons = console["partingLine"]

    # table headline
    headline = "Allen-Dynes-McMillan Formula";
    delimiter = "|"

    # table Results
    ADvalues = [(round(ML_Tc, digits=2)) (round(AD_Tc, digits=2)) (round(BCS_gap, digits=2)) (round(lambda, digits=2)) (round(omega_log, digits=2))]

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
    print(@bold headline)
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
        print(@bold value)
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

    # write everything to inp.log_file
    if inp.flag_log == 1
        print(inp.log_file, logText*"\n")
    end
end


"""
    printSummary()

Summarize Tc calculation and print it to the console
"""
function printSummary(inp, temps, Delta0)

    text = ""
    if all(isnan.(Delta0))
        if minimum(temps) <= 0.5
            text = text*"\n - " * inp.material * " is not a superconductor above T = 0.5 K"
        else
            text = text*"\n - Couldn't find a superconducting gap in the specified area"
            text = text*"\n - Consider searching below T = " * string(minimum(temps)) * " K"
        end
    else
        if inp.TcSearchMode_flag == 0 && ~isnan(last(Delta0))
            text = text*"\n - " * inp.material * " is a superconductor"
            text = text*"\n - Highest given temperature reached"
            text = text*"\n - Tc > " * string(maximum(temps[.~isnan.(Delta0)])) * " K"
        else
            text = text*"\n - " * inp.material * " is a superconductor"
            text = text*"\n - Tc = " * string(maximum(temps[.~isnan.(Delta0)])) * " K!"
        end
    end
    printstyled("\nSummary:", bold=true)
    println(text*"\n")

    if inp.flag_log == 1
        print(inp.log_file, "\nSummary:"*text*"\n\n")
    end

end


"""
    printTextCentered(text, hline[, boldFlag, blanks, delimiter])

print a text centered within a line consisting of delimiters
"""
function printTextCentered(inp, text, hline, boldFlag = false, blanks=3, delimiter = "-")

    lenLeft = Int(ceil((length(hline) - length(text))/2) - blanks)
    lenRight = Int(floor((length(hline) - length(text))/2) - blanks)

    print("\n"*delimiter^lenLeft*" "^blanks)
    if boldFlag
        print(@bold text)
    else
        print(text)
    end
    print(" "^blanks * delimiter^lenRight*"\n")

    if inp.flag_log == 1
        print(inp.log_file, "\n"*"-"^lenLeft*" "^blanks)
        print(inp.log_file, text)
        print(inp.log_file, " "^blanks * "-"^lenRight*"\n")
    end

end


"""
    printTableHeader(console)

Initialize the Table header and print it to the console
"""
function printTableHeader(inp, console)

    width = console["width"]
    header = console["header"]
    initValues = console["InitValues"]

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
        print(@bold value)
        print(delimiter)

        # save for log file
        tableHeader = tableHeader*string(value)*delimiter
    end

    ### parting line ###
    tableHeader = tableHeader*"\n"*replace(tableHline, "." => "|")
    print("\n"*replace(tableHline, "." => "|")*"\n")

    ### Initial values ###
    initValues, strFormat, format = formatTableRow(initValues, width, 2)
    for i in axes(strFormat, 1)
        # print
        if isnothing(initValues[i])
            print(strFormat[i])
        else
            Printf.format(stdout, Printf.Format(strFormat[i]), format[i, 1], " ", format[i, 2], format[i, 3], initValues[i], format[i, 4], " ")
        end
    end

    # write everything to log_file
    if inp.flag_log == 1
        print(inp.log_file,tableHeader*"\n")
        for i in axes(strFormat, 1)
            # print
            if isnothing(initValues[i])
                print(inp.log_file, strFormat[i])
            else
                Printf.format(inp.log_file, Printf.Format(strFormat[i]), format[i, 1], " ", format[i, 2], format[i, 3], initValues[i], format[i, 4], " ")
            end
        end
    end


    console["Hline"] = tableHline

    return console

end


"""
    formatTableHeader(console)

Format the header of the console output s.t. each column has the 
specified length
"""
function formatTableHeader(console::Dict)
    header = console["header"]
    width = console["width"]

    for k in eachindex(header)
        blanks = (width[k]-length(header[k]))/2
        header[k] = " "^Int(floor(blanks))*header[k]*" "^Int(ceil(blanks))
    end


    console["header"] = header

    return console
end

function formatTableHeader(header, width)

    for k in eachindex(header)
        blanks = (width[k]-length(header[k]))/2
        header[k] = " "^Int(floor(blanks))*header[k]*" "^Int(ceil(blanks))
    end

    return header
end


### Format each Row of console output ###
# change to two functions which differ in 
#   - prec::Int = 5
#   - prec::Vector{Int} ?
function formatTableRow(vec, widthCol, prec=5)
    """
    Format the console output s.t. it is aligned to the header
    This is done by introducing a dynamic width and precision via
    "%*s" (requires Julia 1.10)
    The output can be printed via printf   

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
        prec = repeat([prec], length(vec))
    end

    out = Array{String}(undef, length(vec))
    format = zeros(length(vec), 4)
    for k in eachindex(vec)
        value = vec[k]
        width = widthCol[k]

        if isnothing(value) || isnan(value)
            numDig = [0,0]
            spacing = (width - (sum(numDig) + 1)) / 2
            format[k, :] = [floor(spacing), numDig[1], numDig[2], ceil(spacing)]
            out[k] = " "^Int(format[k,1])*"-"*" "^Int(format[k,4])*"|"
        else
            value = round(value, digits=prec[k])
            numDig = numDigits(value)
            if (numDig[1] + prec[k]) > width
                numDig[2] = minimum([0, width - numDig[1]])
            else
                numDig[2] = prec[k]
            end

            if k == 1   # nr iteration, no comma
                out[k] = "|%*s%*.*f%*s|"
                spacing = (width - (numDig[1])) / 2
                format[k, :] = [floor(spacing), numDig[1], 0, ceil(spacing)]
            else
                out[k] = "%*s%*.*f%*s|"
                spacing = (width - (sum(numDig) + 1)) / 2
                format[k, :] = [floor(spacing), numDig[1], numDig[2], ceil(spacing)]
            end
        end
    end

    out[end] = out[end]*"\n"

    return vec, out, format        

end


"""
    numDigits(x)

Gives the number of digits before and after the comma
"""
function numDigits(x)
    digits = zeros(length(x), 2)
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
function printFlagsAsText(inp)
    text = ""
    if inp.material != "Material"
        text *=  " - Material: "*inp.material*" \n"
    end
    # search mode
    if inp.TcSearchMode_flag == 0
        if length(inp.temps) == 1
            text *= " - Tc search area: "*string(inp.temps[1])*" K\n"
        else
            text *= " - Tc search area: "*string(minimum(inp.temps))*" - "*string(maximum(inp.temps))*" K\n"
        end
    elseif inp.TcSearchMode_flag == 1
        text *= " - Tc search mode activated\n"
    end
    
    # cut off
    text *= " - Matsubara cutoff: "*string(inp.omega_c)*" meV\n"

    # fsthick
    if inp.cDOS_flag == 0 
        if inp.fsthick == -1
            text *= " - fsthick: full dos\n"
        else
            text *= " - fsthick: "*string(inp.fsthick)*" meV\n"
        end
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
        text *= " - Morel-Anderson pseudopotential, μ*_AD = "*string(round(inp.muc_AD, digits=3))*" , μ*_ME = "*string(round(inp.muc_ME, digits=3))*"\n"
    end

    print(text)
    if inp.flag_log == 1
        print(inp.log_file, text)
    end

end


"""
    writeToOutFile(Tc ,inp, out_vars, header)

Save results of each iteration and input parameters in a file Summary.txt
"""
function writeToOutFile(Tc ,inp, out_vars, header)
    # write to output file
    name = ""
    if inp.material != "Material"
        name = name*inp.material*"_"
    end
    name = name*"Summary.txt"

    if isfile(inp.outdir*name)
        rm(inp.outdir*name)
    end
    outfile = open(inp.outdir*name, "w")
    
    ### write Tc
    out = "Tc = "*string(Tc)*"\n\n"

    ### calc width of each column
    width = minimum([length.(header) .+ 2])
    width[width .< 12] .= 12
    header = formatTableHeader(header, width)

    ### Define boundary ###
    tableHline = ""
    for w in width
        tableHline = tableHline*"."*"-"^w
    end
    tableHline = tableHline*"."

    ### Define table header ###
    delimiter = "|"
    out = out*tableHline*"\n"*delimiter
    for k in eachindex(header)
        value = header[k]

        out = out*string(value)*delimiter
    end

    ### parting line ###
    out = out*"\n"*replace(tableHline, "." => "|")*"\n"

    ### write to outfile ###
    print(outfile, out)

    ### write values to out file ###
    for i in axes(out_vars, 1)
        var = out_vars[i,:]
        var, strFormat, format = formatTableRow(var, width, 2)
        for j in eachindex(var)
            # print
            if isnan(var[j])
                print(outfile, strFormat[j])
            else
                Printf.format(outfile, Printf.Format(strFormat[j]), format[j, 1], " ", format[j, 2], format[j, 3], var[j], format[j, 4], " ")
            end
        end
    end

    print(outfile, replace(tableHline, "." => " ")*"\n")

    ### Input parameters ###
    print(outfile, "\n\n"*join(inp.all, "\n"))

    close(outfile)
end



"""
    saveSelfEnergyComponents(inp, iwn, Delta, Z; epsilon=nothing,  chi=nothing, phiph=nothing, phic=nothing)

Save the self-energy components into separate files
"""
function saveSelfEnergyComponents(itemp, inp, iwn, Delta, Z; epsilon=nothing,  chi=nothing, phiph=nothing, phic=nothing)

    folder = inp.outdir*"SelfEnergy/"

    if ~isdir(folder)
        mkdir(folder)
    end

    # Z
    open(folder*"Z_"*string(itemp)*"K.dat", "w") do io
        write(io, "#  iωₙ  Z(iωₙ) \n")
        writedlm(io, [iwn Z])
    end

    # Delta 
    if inp.include_Weep == 1
        open(folder*"Delta_"*string(itemp)*"K.dat", "w") do io
            write(io, "# Δ(ϵ, iωₙ) \n")
            writedlm(io, Delta)
        end
    elseif inp.include_Weep == 0
        open(folder*"Delta_"*string(itemp)*"K.dat", "w") do io
            write(io, "#  iωₙ  Δ(iωₙ) \n")
            writedlm(io, [iwn Delta])
        end
     end

    # Chi
    if ~isnothing(chi)
        open(folder*"Chi_"*string(itemp)*"K.dat", "w") do io
            write(io, "#  iωₙ  Χ(iωₙ) \n")
            writedlm(io, [iwn chi])
        end
    end

    # Phiph
    if ~isnothing(phiph)
        open(folder*"Phiph_"*string(itemp)*"K.dat", "w") do io
            write(io, "#  iωₙ  Φp(iωₙ) \n")
            writedlm(io, [iwn phiph])
        end
    end

    # Phic
    if ~isnothing(phic)
        open(folder*"Phic_"*string(itemp)*"K.dat", "w") do io
            write(io, "#  ϵ  Φc(ϵ) \n")
            writedlm(io, [epsilon phic])
        end
    end
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