"""
    Example file

    Several cases are provided as example:
        - 1: cDOS + μ* Tc search - minimal example
        - 2: cDOS + μ* Tc search - different smearing column and typical electronic energy
        - 3: cDOS + μ* explicit temperatures - different smearing column and typical electronic energy
        - 4: vDOS + μ* Tc search - different smearing column and typical electronic energy
        - 5: cDOS + W  Tc search - different smearing column
        - 6: vDOS + W  Tc search - different smearing column

    Cases 7-9 are solved directly on the real frequency axis (RealAxisSolver):
        - 7: cDOS + μ* at a given temperature
        - 8: vDOS + μ* at a given temperature
        - 9: vDOS + W  at a given temperature

"""

using IsoME


## Inputs
# output directory, we recommend to change it
outdir = joinpath(@__DIR__, "output")

# smearing column
ind_smear = 15

# typical electronic energy, used to calculate μ* from μ
typEl = 10000


## Cases
case = 1


if case == 1
    inp = arguments(
        a2f_file    = joinpath(@__DIR__, "Nb.a2f"),
        outdir      = outdir,
    )

elseif case == 2
    inp = arguments(
        a2f_file    = joinpath(@__DIR__, "Nb.a2f"),
        outdir      = outdir,
        ind_smear   = ind_smear,
        typEl       = typEl,
    )

elseif case == 3
    inp = arguments(
        a2f_file    = joinpath(@__DIR__, "Nb.a2f"),
        outdir      = outdir,
        ind_smear   = ind_smear,
        typEl       = typEl,
        temps       = collect(4:2:20)
    )

elseif case == 4
    inp = arguments(
        a2f_file    = joinpath(@__DIR__, "Nb.a2f"),
        dos_file    = joinpath(@__DIR__, "Nb.dos"),
        outdir      = outdir,
        ind_smear   = ind_smear,
        typEl       = typEl,
        cDOS_flag   = 0,
    )

elseif case == 5
    inp = arguments(
        a2f_file        = joinpath(@__DIR__, "Nb.a2f"),
        dos_file        = joinpath(@__DIR__, "Nb.dos"),
        Weep_file       = joinpath(@__DIR__, "Weep.dat"),
        outdir          = outdir,
        ind_smear       = ind_smear,
        include_Weep    = 1,
        cDOS_flag       = 1,
    )

elseif case == 6
    inp = arguments(
        a2f_file        = joinpath(@__DIR__, "Nb.a2f"),
        dos_file        = joinpath(@__DIR__, "Nb.dos"),
        Weep_file       = joinpath(@__DIR__, "Weep.dat"),
        outdir          = outdir,
        ind_smear       = ind_smear,
        include_Weep    = 1,
        cDOS_flag       = 0,
    )

## ----- real axis ----- ##
# The real-axis solver is far more sensitive to its grids than the imaginary-axis one:
# omega_c bounds the ω- and ω'-grids, domega is their step, and depsilon the step of the
# ε-grid in the vDOS cases. See the Real Axis Solver page of the documentation.

elseif case == 7
    inp = arguments(
        a2f_file    = joinpath(@__DIR__, "Nb.a2f"),
        outdir      = outdir,
        ind_smear   = ind_smear,
        typEl       = typEl,
        cDOS_flag   = 1,
        temps       = [6],
        omega_c     = 2000.0,
        domega      = 0.4,
        mixing_beta = 0.5,
        mu          = 0.38,
        conv_thr    = 1e-3,
    )

elseif case == 8
    inp = arguments(
        a2f_file    = joinpath(@__DIR__, "Nb.a2f"),
        dos_file    = joinpath(@__DIR__, "Nb.dos"),
        outdir      = outdir,
        ind_smear   = ind_smear,
        typEl       = typEl,
        cDOS_flag   = 0,
        temps       = [6],
        omega_c     = 2000.0,
        domega      = 0.4,
        depsilon    = 10,
        mixing_beta = 0.5,
        mu          = 0.38,
        conv_thr    = 1e-3,
    )

elseif case == 9
    inp = arguments(
        a2f_file        = joinpath(@__DIR__, "Nb.a2f"),
        dos_file        = joinpath(@__DIR__, "Nb.dos"),
        Weep_file       = joinpath(@__DIR__, "Weep.dat"),
        outdir          = outdir,
        ind_smear       = ind_smear,
        typEl           = typEl,
        cDOS_flag       = 0,
        include_Weep    = 1,
        temps           = [6],
        omega_c         = 2000.0,
        domega          = 0.4,
        depsilon        = 10,
        mixing_beta     = 0.5,
        conv_thr        = 1e-3,
    )
end


## Start the solver
# cases 1-6 are solved on the imaginary axis, cases 7-9 directly on the real axis
if case <= 6
    EliashbergSolver(inp)
else
    RealAxisSolver(inp)
end
