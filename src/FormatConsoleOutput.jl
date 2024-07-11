"""
Format the console output nicely
    - Start message
    - Iteration output

Julia Packages:
    - 

Comments:
    - 

"""


function printStartMessage(console)
    """
    Start message

    -------------------------------------------------------------------
    Input:
        -

    --------------------------------------------------------------------
    Output:
        -

    --------------------------------------------------------------------
    Comments:
        - 

    --------------------------------------------------------------------
    """

    strIsoMe = "\n\n"
    strIsoMe = strIsoMe*"   _                 __  __   ______ \n"
    strIsoMe = strIsoMe*"  (_)               |  \\/  | |  ____|\n"
    strIsoMe = strIsoMe*"   _   ___    ___   | \\  / | | |__   \n"
    strIsoMe = strIsoMe*"  | | / __|  / _ \\  | |\\/| | |  __|  \n"
    strIsoMe = strIsoMe*"  | | \\__ \\ | (_) | | |  | | | |____ \n"
    strIsoMe = strIsoMe*"  |_| |___/  \\___/  |_|  |_| |______|\n\n"

    #println("\nAuthors: ")
    #println("  - Christoph Heil")
    #println("  - Eva Kogler")
    #println("  - Dominik Spath\n")
    #println("Version: 1.0")

    strLine = "-"^(sum(console["width"])+length(console["width"])+1)
    strEliash = "Eliashberg Solver started"

    print(strIsoMe)
    print(strLine)
    if flag_log == 1
        print(log_file, strIsoMe)
        print(log_file, strLine)
    end
    printTextCentered(strEliash, strLine, true)
    print(strLine*"\n\n\n")
    if flag_log == 1
        print(log_file, strLine*"\n\n\n")
    end

    console["partingLine"] = strLine
    
    return console
                                    
end

### 
function printADtable(console)
    """
    Start message

    -------------------------------------------------------------------
    Input:
        - partingLineCons: horizontal line console

    --------------------------------------------------------------------
    Output:
        -

    --------------------------------------------------------------------
    Comments:
        - 

    --------------------------------------------------------------------
    """

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
    logText = logText*" "^blanksAD*replace(Hline, "." => ":")*"\n"
    println(" "^blanksAD*replace(Hline, "." => ":"))

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

    # write everything to log_file
    if flag_log == 1
        print(log_file, logText*"\n")
    end
end


### Summary of complete calculation ###
function printSummary()
    """
    Start message

    -------------------------------------------------------------------
    Input:
        -

    --------------------------------------------------------------------
    Output:
        -

    --------------------------------------------------------------------
    Comments:
        - 

    --------------------------------------------------------------------
    """


    text = ""
    if all(isnan.(Delta0))
        if minimum(temps) <= 1
            text = text*"\n - " * material * " is not a superconductor"
        else
            text = text*"\n - Couldn't find a superconducting gap in the specified area"
            text = text*"\n - Consider searching below T = " * string(minimum(temps)) * " K"
        end
    else
        if TcSearchMode_flag == 0 && ~isnan(last(Delta0))
            text = text*"\n - " * material * " is a superconductor"
            text = text*"\n - Highest given temperature reached"
            text = text*"\n - Tc > " * string(temps[findlast(.~isnan.(Delta0))]) * " K"
        else
            text = text*"\n - " * material * " is a superconductor"
            text = text*"\n - Tc = " * string(temps[findlast(.~isnan.(Delta0))]) * " K!"
        end
    end
    printstyled("\nSummary:", bold=true)
    println(text*"\n")

    if flag_log == 1
        print(log_file, "\nSummary:"*text*"\n\n")
    end

end

### Print a text centered within parting line ###
function printTextCentered(text, hline, boldFlag = false, blanks=3, del = "-")
    """
    print a text centered within a parting line

    -------------------------------------------------------------------
    Input:
        - text:     text
        - hline:    horizontal line 
        - boldFlag: text should be bold
        - blanks:   blanks between hline and text
        - del:      delimiter that should be used to fill up the space

    --------------------------------------------------------------------
    Output:
        - 

    --------------------------------------------------------------------
    Comments:
        - 

    --------------------------------------------------------------------
    """


    lenLeft = Int(ceil((length(hline) - length(text))/2) - blanks)
    lenRight = Int(floor((length(hline) - length(text))/2) - blanks)

    print("\n"*del^lenLeft*" "^blanks)
    if boldFlag
        print(@bold text)
    else
        print(text)
    end
    print(" "^blanks * del^lenRight*"\n")

    if flag_log == 1
        print(log_file, "\n"*"-"^lenLeft*" "^blanks)
        print(log_file, text)
        print(log_file, " "^blanks * "-"^lenRight*"\n")
    end

end


### print header of table to console ###
function printTableHeader(console)
    """
    Initialize the Table header for the console output and print it
    Return it for log file 

    -------------------------------------------------------------------
    Input:
        - header:       raw header                          | string
        - initValues:   values of initial iteration         | vector
        - width:  width of each column in the header  | vector

    --------------------------------------------------------------------
    Output:
        - str:      formated table header       | string

    --------------------------------------------------------------------
    Comments:
        -

    --------------------------------------------------------------------
    """


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
    if flag_log == 1
        print(log_file,tableHeader*"\n")
        for i in axes(strFormat, 1)
            # print
            if isnothing(initValues[i])
                print(log_file, strFormat[i])
            else
                Printf.format(log_file, Printf.Format(strFormat[i]), format[i, 1], " ", format[i, 2], format[i, 3], initValues[i], format[i, 4], " ")
            end
        end
    end


    console["Hline"] = tableHline

    return console

end


### Format header of console output ###
function formatTableHeader(console::Dict)
    """
    Format the header of the table s.t. each column has the specified
    length

    -------------------------------------------------------------------
    Input:
        - header:    raw header                          | string
        - width:     width of each column in the header  | vector

    --------------------------------------------------------------------
    Output:
        - str:      formated table header       | string

    --------------------------------------------------------------------
    Comments:
        -

    --------------------------------------------------------------------
    """

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
    """
    Format the header of the table s.t. each column has the specified
    length

    -------------------------------------------------------------------
    Input:
        - header:    raw header                          | string
        - width:     width of each column in the header  | vector

    --------------------------------------------------------------------
    Output:
        - str:      formated table header       | string

    --------------------------------------------------------------------
    Comments:
        -

    --------------------------------------------------------------------
    """

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

        if isnothing(value)
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


function numDigits(x)
    """
    Gives the number of digits before and after the comma

    -------------------------------------------------------------------
    Input:
        x:      vector

    --------------------------------------------------------------------
    Output:
        digits:  digits before/after comma for each element in x

    --------------------------------------------------------------------
    Comments:
        - 
    --------------------------------------------------------------------
    """

    digits = zeros(length(x), 2)
    for k in eachindex(x)
        str = split(string(x[k]), ".")
        digits[k,:] = [length(str[1]), length(str[2])]
    end
    return digits

end


# Format flags to text for console
function printFlagsAsText()
    """
    Print the flag values as text to the console

    -------------------------------------------------------------------
    Input:

    --------------------------------------------------------------------
    Output:
        text:  console output

    --------------------------------------------------------------------
    Comments:
        - 
    --------------------------------------------------------------------
    """

    text =  " - Material: "*material*" \n"
    # search mode
    if TcSearchMode_flag == 0
        text *= " - Tc search area: "*string(minimum(temps))*" - "*string(maximum(temps))*" K\n"
    elseif TcSearchMode_flag == 1
        text *= " - Tc auto search\n"
    end
    
    # cut off
    text *= " - Matsubara cutoff: "*string(omega_c)*" meV\n"
    
    # cDos
    if cDOS_flag == 0
        if mu_flag == 1
            text *= " - Variable DoS with μ-update\n"
        else
            text *= " - Variable DoS with constant μ = ϵ_F\n"
        end
    elseif cDOS_flag == 1
        text *= " - Constant DoS\n"
    end

    # Weep
    if include_Weep == 1
        text *= " - Full Coulomb interaction given in "*unitWeepFile*"\n"
    elseif include_Weep == 0
        text *= " - Anderson pseudopotential, μ* = "*string(muc)*" , μ*_ME = "*string(round(muc_ME, digits=2))*"\n"
    end

    # other DoS
    #if other_dos_flag == 1
    #    text *= " - Other Dos given in "*unitScfDosFile*"\n"
    #else
    #    text *= " - Dos given in "*unitDosFile*"\n"
    #end
    #text = text*"\n"


    print(text)
    if flag_log == 1
        print(log_file, text)
    end

end

