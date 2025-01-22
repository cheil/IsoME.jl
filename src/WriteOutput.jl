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

Start message
"""
function printStartMessage(console, log_file)

    strLine = "-"^(sum(console["width"])+length(console["width"])+1)

    strAuthors = "  Authors: Christoph Heil, Eva Kogler, Dominik Spath\n\n"

    strEliash = "Eliashberg Solver started"

    #printTextCentered(log_file, "Version 1.0", strLine, false)
    print(strAuthors)
    print(strLine)
    # log file
    print(log_file, strAuthors)
    print(log_file, strLine)

    printTextCentered(strEliash, strLine, file = log_file, bold = true)
    print(strLine*"\n\n\n")
    # log file
    print(log_file, strLine*"\n\n\n")
    
    console["partingLine"] = strLine
    
    return console
                                    
end


"""
    printADtable(console)

Print the Allen-Dynes results to the console.
"""
function printADtable(console, ML_Tc, AD_Tc, BCS_gap, lambda, omega_log, log_file)

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
    println(text)

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
            print(@bold text)
        else
            print(text)
        end
        print(rightText)
    end

    # print to file
    if isfile(file)
        print(file, leftText)
        print(file, text)
        print(file, rightText)
    end

end


"""
    printTableHeader(console)

Initialize the Table header and print it to the console
"""
function printTableHeader(console, log_file)

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
    print(log_file, tableHeader * "\n")
    for i in axes(strFormat, 1)
        # print
        if isnothing(initValues[i])
            print(log_file, strFormat[i])
        else
            Printf.format(log_file, Printf.Format(strFormat[i]), format[i, 1], " ", format[i, 2], format[i, 3], initValues[i], format[i, 4], " ")
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
function formatTableRow(vec, widthCol, prec=5, logConsole=true)
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

            if k == 1    # nr iteration, no comma
                out[k] = "|%*s%*.*f%*s|"
                if logConsole
                    spacing = (width - (numDig[1])) / 2
                    format[k, :] = [floor(spacing), numDig[1], 0, ceil(spacing)]
                else
                    spacing = (width - (sum(numDig) + 1)) / 2
                    format[k, :] = [floor(spacing), numDig[1], numDig[2], ceil(spacing)]
                end
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
function printFlagsAsText(inp, log_file)
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
    text *= " - Matsubara cutoff: "*string(inp.omega_c)*" meV\n"

    # encut
    if inp.cDOS_flag == 0 
        if inp.encut == -1
            text *= " - encut: full dos\n"
        else
            text *= " - encut: "*string(inp.encut)*" meV\n"
        end
    end

    # shiftcut
    if inp.cDOS_flag == 0 
        if inp.shiftcut == -1
            text *= " - shiftcut: full dos\n"
        else
            text *= " - shiftcut: "*string(inp.shiftcut)*" meV\n"
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
        text *= " - Morel-Anderson pseudopotential\n"
        text *= "     - μ*_AD = "*string(round(inp.muc_AD, digits=3))*"\n"
        text *= "     - μ*_ME = "*string(round(inp.muc_ME, digits=3))*"\n"
    end

    print(text)
    
    # log file
    print(log_file, text)


end


"""
"""
function printTee(log_file, text)
    print(text)
    print(log_file, text)
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
    flush(log_file)
    close(log_file) 

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
    writeToCrashFile(inp)

Save exception in CRASH file
"""
function writeToCrashFile(inp)
    crashFile = open(inp.outdir * "CRASH", "a")
    print(crashFile, current_exceptions())
    print(crashFile, "\n\n")
    close(crashFile)
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
    print(outfile, replace(replace(replace(join(inp.all, "\n"), "-1.0"=>"-"), "Number[-1]"=>"-"), "-1"=>"-"))

    close(outfile)
end


"""
    summarizeResults(Tc, out_vars, header)

Save Delta(0) at each temperature
"""
function createSummaryFile(inp, Tc, out_vars, header)
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
      
      # save header & gap
      #writedlm(outfile, header)
      writedlm(outfile, round.(out_vars, digits=2), '\t')

      close(outfile)

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
    xlim_max = round(maximum(a2f_omega_fine) / 10 * 1.01, RoundUp) * 10
    xtick_val = 0:10:xlim_max
    ylim_max = round(maximum(a2f_fine), RoundUp)

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
        print(@blue "Info: ")
        println("No superconducting gap found - skipping plot\n")

        print(log_file,  "Info: No superconducting gap found - skipping plot\n")
    else
        # print gap vs. temperature
        Delta0_plot = Delta0[.~isnan.(Delta0)]
        temps_plot = temps[.~isnan.(Delta0)]
        order = sortperm(temps_plot)
        temps_plot = temps_plot[order]
        Delta0_plot = Delta0_plot[order]

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
        
        # normalize x 
        color_values = (temps_plot)/max_gradient_val/2
        
        # Map colors: Use blue gradient for values <= max_gradient_val, red gradient for rest
        marker_colors = [v < 1.0 ? my_gradient[v] : my_gradient[v-1] for v in color_values]
        
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
                if all(DeltaSecDer(temps_fit).<0) && (TcFit > Tc[1]) #&& (TcFit < Tc[2] || isnan(Tc[2])) && TcFit < 1e3 # plot only if TcFit is reasonable ??
                    plot!(temps_fit, Delta_fit, linestyle=:dash, linewidth=1, color=:gray, z_order=:back)
            
                    # adjust limits
                    xlim_max = maximum([xlim_max, TcFit+1])
                    ylim_max = maximum([ylim_max, Delta0Fit+1])
            
                    if xlim_max <= 10
                        xtick_val = 0:1:xlim_max
                    elseif xlim_max > 1e3
                        error("TcFit unreasonbale")
                    elseif xlim_max <= 20
                        xtick_val = 0:2:xlim_max
                    else
                        xtick_val = 0:10:xlim_max
                    end
                end
            catch  ex
                #rethrow(ex)
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
        write(io, "#  iωₙ / meV      Z(iωₙ) / 1 \n")
        writedlm(io, [iwn Z], '\t')
    end

    # Delta 
    if inp.include_Weep == 1
        open(folder*"Delta_"*string(itemp)*"K.dat", "w") do io
            write(io, "# ε / meV      iωₙ / meV      Δ(ϵ, iωₙ) / meV \n")
            writedlm(io, [repeat(epsilon, inner=length(iwn)) repeat(iwn, length(epsilon)) Delta'[:]], '\t')
        end
        
        # w/o epsilon and mat.freqs.
        #open(folder*"Delta_"*string(itemp)*"K.dat", "w") do io
        #    write(io, "# Δ(ϵ, iωₙ) / meV \n")
        #    writedlm(io, [Delta], '\t')
        #end
    elseif inp.include_Weep == 0
        open(folder*"Delta_"*string(itemp)*"K.dat", "w") do io
            write(io, "#  iωₙ / meV      Δ(iωₙ) / meV \n")
            writedlm(io, [iwn Delta], '\t')
        end
     end

    # Chi
    if ~isnothing(chi)
        open(folder*"Chi_"*string(itemp)*"K.dat", "w") do io
            write(io, "#  iωₙ / meV      Χ(iωₙ) / meV \n")
            writedlm(io, [iwn chi], '\t')
        end
    end

    # Phiph
    if ~isnothing(phiph)
        open(folder*"Phiph_"*string(itemp)*"K.dat", "w") do io
            write(io, "#  iωₙ /meV      Φp(iωₙ) / meV \n")
            writedlm(io, [iwn phiph], '\t')
        end
    end

    # Phic
    if ~isnothing(phic)
        open(folder*"Phic_"*string(itemp)*"K.dat", "w") do io
            write(io, "#  ϵ / meV      Φc(ϵ) / meV \n")
            writedlm(io, [epsilon phic], '\t')
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