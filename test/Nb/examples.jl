# Example file demonstrating how each of the four different modes can be started

using IsoME

case = 5
# 1: cDOS + μ* Tc search
# 2: cDOS + μ* explicit temperatures
# 3: vDOS + μ* Tc search
# 4: cDOS + W Tc search
# 5: vDOS + W Tc search


if case == 1
    inp = arguments(
        a2f_file    = joinpath(@__DIR__, "Nb.a2F"),
        outdir      = joinpath(@__DIR__, "output"),
    )

elseif case == 2
    inp = arguments(
        a2f_file    = joinpath(@__DIR__, "Nb.a2F"),
        outdir      = joinpath(@__DIR__, "output"),
        TcSearchMode_flag = 0,
        temps       = collect(2:2:20) 
    )

elseif case == 3
    inp = arguments(
        a2f_file    = joinpath(@__DIR__, "Nb.a2F"),
        dos_file    = joinpath(@__DIR__, "Nb.dos"),
        outdir      = joinpath(@__DIR__, "output"),
        cDOS_flag   = 0,
    )

elseif case == 4
    inp = arguments(
        a2f_file        = joinpath(@__DIR__, "Nb.a2F"),
        dos_file        = joinpath(@__DIR__, "Nb.dos"),
        Weep_file       = joinpath(@__DIR__, "Weep.dat"),
        Wen_file        = joinpath(@__DIR__, "Wen.dat"),
        outdir          = joinpath(@__DIR__, "output"),
        include_Weep    = 1,
        cDOS_flag       = 1,
    )

elseif case == 5
    inp = arguments(
        a2f_file    = joinpath(@__DIR__, "Nb.a2F"),
        dos_file    = joinpath(@__DIR__, "Nb.dos"),
        Weep_file   = joinpath(@__DIR__, "Weep.dat"),
        Wen_file    = joinpath(@__DIR__, "Wen.dat"),
        outdir      = joinpath(@__DIR__, "output"),
        include_Weep    = 1,
        cDOS_flag       = 0,
    )
end

# Run the solver with the specified input
EliashbergSolver(inp)




