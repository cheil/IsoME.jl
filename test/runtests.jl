# !!! Automatic test !!!
# Add further tests for sparse sampling, encut, interpoaltion, ....

using IsoME
using Test


@testset "IsoME.jl" begin

    ### 1.TEST: Nb cDOS mu* ###
    inp = arguments(
                    a2f_file = joinpath(@__DIR__, "Nb/Nb.a2F"),
                    dos_file = joinpath(@__DIR__, "Nb/Nb.dos"),
                    Weep_file = joinpath(@__DIR__, "Nb/Weep.dat"),
                    flag_figure = 0,
                    returnTc    = true,
                    testMode    = true,
                    outdir      = "./test/Nb/output/",
                    ind_smear   = 15,
                    typEl       = 10000,
                    )

    Tc = EliashbergSolver(inp)
    @test Tc == [9,10]


    ### 2.TEST: Nb vDOS mu* ###
    inp = arguments(
                    a2f_file = joinpath(@__DIR__, "Nb/Nb.a2F"),
                    dos_file = joinpath(@__DIR__, "Nb/Nb.dos"),
                    Weep_file = joinpath(@__DIR__, "Nb/Weep.dat"),
                    flag_figure = 0,
                    returnTc    = true,
                    testMode    = true,
                    outdir      = "./test/Nb/output/",
                    cDOS_flag   = 0,
                    ind_smear   = 15,
                    typEl       = 10000,
                    )

    Tc = EliashbergSolver(inp)
    @test Tc == [8, 9]


    ### 3.TEST: Nb cDOS W ###
    inp = arguments(
        a2f_file    = joinpath(@__DIR__, "Nb/Nb.a2F"),
        dos_file    = joinpath(@__DIR__, "Nb/Nb.dos"),
        Weep_file   = joinpath(@__DIR__, "Nb/Weep.dat"),
        flag_figure =0,
        returnTc    = true,
        testMode    = true,
        outdir      = "./test/Nb/output/",
        cDOS_flag   = 1,
        include_Weep = 1,
        ind_smear   = 15,
        typEl       = 10000,
    )

    Tc = EliashbergSolver(inp)
    @test Tc == [7, 8]


    ### 4.TEST: Nb vDOS W ###
    inp = arguments(
        a2f_file = joinpath(@__DIR__, "Nb/Nb.a2F"),
        dos_file = joinpath(@__DIR__, "Nb/Nb.dos"),
        Weep_file = joinpath(@__DIR__, "Nb/Weep.dat"),
        flag_figure=0,
        returnTc    = true,
        testMode    = true,
        outdir      = "./test/Nb/output/",
        cDOS_flag = 0,
        include_Weep = 1,
        ind_smear   = 15,
        typEl       = 10000,
    )

    Tc = EliashbergSolver(inp)
    @test Tc == [7, 8]

end


