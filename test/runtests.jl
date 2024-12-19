# !!! Automatic test - Do not change !!!

using IsoME
using Test


@testset "IsoME.jl" begin

    outdir  = "test/Nb/output/"

    ### 1.TEST: Nb cDOS mu* ###
    inp = arguments(
                    a2f_file = joinpath(@__DIR__, "Nb/Nb.a2F"),
                    dos_file = joinpath(@__DIR__, "Nb/Nb.dos"),
                    Weep_file = joinpath(@__DIR__, "Nb/Weep.dat"),
                    flag_figure = 0,
                    returnTc    = true,
                    outdir = outdir*"cDOS_mu/",
                    )

    Tc = EliashbergSolver(inp)
    @test Tc == [13,14]


    ### 2.TEST: Nb vDOS mu* ###
    inp = arguments(
                    a2f_file = "test/Nb/Nb.a2F",
                    dos_file = "test/Nb/Nb.dos",
                    Weep_file = "test/Nb/Weep.dat",
                    Wen_file = "test/Nb/Wen.dat",
                    flag_figure = 0,
                    returnTc    = true,
                    outdir = outdir*"vDOS_mu/",
                    cDOS_flag   = 0,
                    )

    Tc = EliashbergSolver(inp)
    @test Tc == [12, 13]


    ### 3.TEST: Nb cDOS W ###
    inp = arguments(
        a2f_file="test/Nb/Nb.a2F",
        dos_file="test/Nb/Nb.dos",
        Weep_file="test/Nb/Weep.dat",
        Wen_file="test/Nb/Wen.dat",
        flag_figure=0,
        returnTc=true,
        outdir=outdir * "cDOS_W/",
        cDOS_flag = 1,
        include_Weep = 1,
    )

    Tc = EliashbergSolver(inp)
    @test Tc == [13, 14]


    ### 4.TEST: Nb vDOS W ###
    inp = arguments(
        a2f_file="test/Nb/Nb.a2F",
        dos_file="test/Nb/Nb.dos",
        Weep_file="test/Nb/Weep.dat",
        Wen_file="test/Nb/Wen.dat",
        flag_figure=0,
        returnTc=true,
        outdir=outdir * "vDOS_W/",
        cDOS_flag = 0,
        include_Weep = 1,
    )

    Tc = EliashbergSolver(inp)
    @test Tc == [13, 14]


    # remove outdir
    rm(outdir, recursive=true)

end


