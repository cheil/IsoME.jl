"""
Example file demonstrating how each of the four different modes can be started
"""
using IsoME

case = 1
# 1: cDOS + μ*
# 2: vDOS + μ*
# 3: cDOS + W
# 4: vDOS + W


if case == 1
    inp = arguments(
        a2f_file    = "./Nb.a2F",
        outdir      = "./output/",
    )

    EliashbergSolver(inp)

elseif case == 2
    inp = arguments(
        a2f_file    = "./Nb.a2F",
        dos_file    = "./Nb.dos",
        outdir      = "./output/",
        cDOS_flag   = 0,
    )

    EliashbergSolver(inp)

elseif case == 3
    inp = arguments(
        a2f_file        = "./Nb.a2F",
        dos_file        = "./Nb.dos",
        Weep_file       = "./Weep.dat",
        Wen_file        = "./Wen.dat",
        outdir          = "./output/",
        include_Weep    = 1,
        cDOS_flag       = 1,
    )

    EliashbergSolver(inp)

elseif case == 4
    inp = arguments(
        a2f_file    = "./Nb.a2F",
        dos_file    = "./Nb.dos",
        Weep_file   = "./Weep.dat",
        Wen_file    = "./Wen.dat",
        outdir      = "./output/",
        include_Weep    = 1,
        cDOS_flag       = 0,
    )

    EliashbergSolver(inp)

end




