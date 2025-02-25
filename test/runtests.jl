# !!! Automatic test - Do not change !!!
# Add further tests for sparse sampling, encut, interpoaltion, ....

using IsoME
using Test


@testset "IsoME.jl" begin

    ### 1.TEST: Nb cDOS mu* ###
    inp = arguments(
                    a2f_file = joinpath(@__DIR__, "Nb/Nb.a2F"),
                    dos_file = joinpath(@__DIR__, "Nb/Nb.dos"),
                    Weep_file = joinpath(@__DIR__, "Nb/Weep.dat"),
                    Wen_file = joinpath(@__DIR__, "Nb/Wen.dat"),
                    flag_figure = 0,
                    returnTc    = true,
                    testMode    = true,
                    outdir = "",
                    muc_AD = 0.12,
                    )

    Tc = EliashbergSolver(inp)
    @test Tc == [13,14]


    ### 2.TEST: Nb vDOS mu* ###
    inp = arguments(
                    a2f_file = joinpath(@__DIR__, "Nb/Nb.a2F"),
                    dos_file = joinpath(@__DIR__, "Nb/Nb.dos"),
                    Weep_file = joinpath(@__DIR__, "Nb/Weep.dat"),
                    Wen_file = joinpath(@__DIR__, "Nb/Wen.dat"),
                    flag_figure = 0,
                    returnTc    = true,
                    testMode    = true,
                    outdir = "",
                    cDOS_flag   = 0,
                    muc_AD = 0.12,
                    )

    Tc = EliashbergSolver(inp)
    @test Tc == [12, 13]


    ### 3.TEST: Nb cDOS W ###
    inp = arguments(
        a2f_file = joinpath(@__DIR__, "Nb/Nb.a2F"),
        dos_file = joinpath(@__DIR__, "Nb/Nb.dos"),
        Weep_file = joinpath(@__DIR__, "Nb/Weep.dat"),
        Wen_file = joinpath(@__DIR__, "Nb/Wen.dat"),
        flag_figure=0,
        returnTc    = true,
        testMode    = true,
        outdir = "",
        cDOS_flag = 1,
        include_Weep = 1,
    )

    Tc = EliashbergSolver(inp)
    @test Tc == [13, 14]


    ### 4.TEST: Nb vDOS W ###
    inp = arguments(
        a2f_file = joinpath(@__DIR__, "Nb/Nb.a2F"),
        dos_file = joinpath(@__DIR__, "Nb/Nb.dos"),
        Weep_file = joinpath(@__DIR__, "Nb/Weep.dat"),
        Wen_file = joinpath(@__DIR__, "Nb/Wen.dat"),
        flag_figure=0,
        returnTc    = true,
        testMode    = true,
        outdir = "",
        cDOS_flag = 0,
        include_Weep = 1,
    )

    Tc = EliashbergSolver(inp)
    @test Tc == [13, 14]

end


