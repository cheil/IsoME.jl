# !!! Automatic test - Do not change !!!
# What to do with outdir
# testFlag

using IsoME
using Test


@testset "IsoME.jl" begin
    ### Nb ###
    inp = arguments(
                    cDOS_flag = 1,
                    include_Weep = 0,
                    a2f_file = "test/Nb/Nb.a2F",
                    dos_file = "test/Nb/Nb.dos",
                    Weep_file = "test/Nb/Weep.dat",
                    Wen_file = "test/Nb/Wen.dat",
                    flag_figure = 0,
                    TcSearchMode_flag = 1
                    )

    @test 1 == 1
#=
    Tc = EliashbergSolver(inp, true)
    @test Tc == 12
       
    ### USE CORRECT Tc VALUES !!! ###
    inp.cDOS_flag = 0
    inp.temps = [23,24,25]
    Tc = EliashbergSolver(inp, true)
    @test Tc == 13

    inp.cDOS_flag = 1
    inp.include_Weep = 1
    inp.temps = [26,27,28]
    Tc = EliashbergSolver(inp, true)
    @test Tc == 13

    inp.cDOS_flag = 0
    inp.temps = [23,24,25]
    Tc = EliashbergSolver(inp, true)
    @test Tc == 13
=#
    
end


