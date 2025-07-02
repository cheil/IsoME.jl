"""
    Example file 

    Several cases are provided as example:
        - 1: cDOS + μ* Tc search - minimal example
        - 2: cDOS + μ* Tc search - different smearing column and typical electronic energy
        - 3: cDOS + μ* explicit temperatures - different smearing column and typical electronic energy
        - 4: vDOS + μ* Tc search - different smearing column and typical electronic energy
        - 5: cDOS + W  Tc search - different smearing column 
        - 6: vDOS + W  Tc search - different smearing column

"""

using IsoME


# Inputs
# output directory, we recommend to change it
outdir = joinpath(@__DIR__, "output")

# smearing column
ind_smear = 15

# typical electronic energy, used to calculate μ* from μ
typEl = 10000


# Cases
case = 1


if case == 1
    inp = arguments(
        a2f_file    = joinpath(@__DIR__, "Nb.a2F"),
        outdir      = outdir,
    )

elseif case == 2
    inp = arguments(
        a2f_file    = joinpath(@__DIR__, "Nb.a2F"),
        outdir      = outdir,
        ind_smear   = ind_smear,
        typEl       = typEl,
    )

elseif case == 3
    inp = arguments(
        a2f_file    = joinpath(@__DIR__, "Nb.a2F"),
        outdir      = outdir,
        ind_smear   = ind_smear,
        typEl       = typEl,
        temps       = collect(4:2:20) 
    )

elseif case == 4
    inp = arguments(
        a2f_file    = joinpath(@__DIR__, "Nb.a2F"),
        dos_file    = joinpath(@__DIR__, "Nb.dos"),
        outdir      = outdir,
        ind_smear   = ind_smear,
        typEl       = typEl,
        cDOS_flag   = 0,
    )

elseif case == 5
    inp = arguments(
        a2f_file        = joinpath(@__DIR__, "Nb.a2F"),
        dos_file        = joinpath(@__DIR__, "Nb.dos"),
        Weep_file       = joinpath(@__DIR__, "Weep.dat"),
        outdir          = outdir,
        ind_smear       = ind_smear,
        include_Weep    = 1,
        cDOS_flag       = 1,
    )

elseif case == 6
    inp = arguments(
        a2f_file        = joinpath(@__DIR__, "Nb.a2F"),
        dos_file        = joinpath(@__DIR__, "Nb.dos"),
        Weep_file       = joinpath(@__DIR__, "Weep.dat"),
        outdir          = outdir,
        ind_smear       = ind_smear,
        include_Weep    = 1,
        cDOS_flag       = 0,
    )
end


# Start the Eliashberg Solver
EliashbergSolver(inp)




