using IsoME
using Test


@testset "IsoME.jl" begin
    ### Nb ###
    inp = arguments(
                    temps = [10, 14],
                    cDOS_flag = 1,
                    include_Weep = 0,
                    a2f_file = "Nb/a2f_k48.dat",
                    dos_file = "Nb/DOS_scf.dat",
                    Weep_file = "Nb/Weep.dat",
                    dosW_file = "Nb/DOSW.dat",
                    flag_log = 0,
                    flag_figure = 0,
                    flag_outfile = 0
                    )

    print(pwd())

    Delta = EliashbergSolver(inp, true)
    @test Delta == [3.14, 2.61]

    #=
    inp.cDOS_flag = 0
    Delta = EliashbergSolver(inp, true)
    @test Delta == [2.31, ]

    inp.cDOS_flag = 1
    inp.include_Weep = 1
    Delta = EliashbergSolver(inp, true)
    @test Delta == [3.14, 2.61]

    inp.cDOS_flag = 0
    Delta = EliashbergSolver(inp, true)
    @test Delta == [3.14, 2.61]
    =# 
end


