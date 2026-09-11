# !!! Automatic test !!!
# Add further tests for sparse sampling, encut, interpoaltion, ....

using IsoME
using Test


@testset "IsoME.jl" begin

    ############################################################
    # ------------------- imaginary axis --------------------- #
    ############################################################

    ### 1.TEST: Nb cDOS mu* ###
    inp = arguments(
                    a2f_file = joinpath(@__DIR__, "Nb/Nb.a2f"),
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


    ### 2.TEST: Nb vDOS mu*, fixed mu ###
    inp = arguments(
                    a2f_file = joinpath(@__DIR__, "Nb/Nb.a2f"),
                    dos_file = joinpath(@__DIR__, "Nb/Nb.dos"),
                    Weep_file = joinpath(@__DIR__, "Nb/Weep.dat"),
                    flag_figure = 0,
                    returnTc    = true,
                    testMode    = true,
                    outdir      = "./test/Nb/output/",
                    cDOS_flag   = 0,
                    mu_flag     = 0,
                    ind_smear   = 15,
                    typEl       = 10000,
                    )

    Tc = EliashbergSolver(inp)
    @test Tc == [7,8]


    ### 3.TEST: Nb vDOS mu*, mu-update ###
    inp = arguments(
                    a2f_file = joinpath(@__DIR__, "Nb/Nb.a2f"),
                    dos_file = joinpath(@__DIR__, "Nb/Nb.dos"),
                    Weep_file = joinpath(@__DIR__, "Nb/Weep.dat"),
                    flag_figure = 0,
                    returnTc    = true,
                    testMode    = true,
                    outdir      = "./test/Nb/output/",
                    cDOS_flag   = 0,
                    mu_flag     = 1,
                    ind_smear   = 15,
                    typEl       = 10000,
                    )

    Tc = EliashbergSolver(inp)
    @test Tc == [7,8]


    ### 4.TEST: Nb cDOS W ###
    # inp = arguments(
    #     a2f_file    = joinpath(@__DIR__, "Nb/Nb.a2f"),
    #     dos_file    = joinpath(@__DIR__, "Nb/Nb.dos"),
    #     Weep_file   = joinpath(@__DIR__, "Nb/Weep.dat"),
    #     flag_figure =0,
    #     returnTc    = true,
    #     testMode    = true,
    #     outdir      = "./test/Nb/output/",
    #     cDOS_flag   = 1,
    #     include_Weep = 1,
    #     ind_smear   = 15,
    #     typEl       = 10000,
    # )

    # Tc = EliashbergSolver(inp)
    # @test Tc == [7, 8]


    ### 5.TEST: Nb vDOS W, fixed mu ###
    inp = arguments(
        a2f_file = joinpath(@__DIR__, "Nb/Nb.a2f"),
        dos_file = joinpath(@__DIR__, "Nb/Nb.dos"),
        Weep_file = joinpath(@__DIR__, "Nb/Weep.dat"),
        flag_figure=0,
        returnTc    = true,
        testMode    = true,
        outdir      = "./test/Nb/output/",
        cDOS_flag = 0,
        include_Weep = 1,
        mu_flag     = 0,
        ind_smear   = 15,
        typEl       = 10000,
    )

    Tc = EliashbergSolver(inp)
    @test Tc == [7, 8]


    ### 6.TEST: Nb vDOS W, mu-update ###
    inp = arguments(
        a2f_file = joinpath(@__DIR__, "Nb/Nb.a2f"),
        dos_file = joinpath(@__DIR__, "Nb/Nb.dos"),
        Weep_file = joinpath(@__DIR__, "Nb/Weep.dat"),
        flag_figure=0,
        returnTc    = true,
        testMode    = true,
        outdir      = "./test/Nb/output/",
        cDOS_flag = 0,
        include_Weep = 1,
        mu_flag     = 1,
        ind_smear   = 15,
        typEl       = 10000,
    )

    Tc = EliashbergSolver(inp)
    @test Tc == [7, 8]


    ############################################################
    # --------------------- real axis ------------------------ #
    ############################################################

    ### 7.TEST: Nb cDOS mu* - real axis ###
    inp = arguments(
        a2f_file    = joinpath(@__DIR__, "Nb/Nb.a2f"),
        dos_file    = joinpath(@__DIR__, "Nb/Nb.dos"),
        Weep_file   = joinpath(@__DIR__, "Nb/Weep.dat"),
        flag_figure = 0,
        returnTc    = true,
        testMode    = true,
        outdir      = "./test/Nb/output/",
        cDOS_flag   = 1,
        ind_smear   = 15,
        typEl       = 10000,
    )

    Tc = RealAxisSolver(inp)
    @test Tc == [9, 10]


    ### 8.TEST: Nb vDOS mu* - real axis, fixed mu ###
    inp = arguments(
        a2f_file    = joinpath(@__DIR__, "Nb/Nb.a2f"),
        dos_file    = joinpath(@__DIR__, "Nb/Nb.dos"),
        Weep_file   = joinpath(@__DIR__, "Nb/Weep.dat"),
        flag_figure = 0,
        returnTc    = true,
        testMode    = true,
        outdir      = "./test/Nb/output/",
        cDOS_flag   = 0,
        mu_flag     = 0,
        ind_smear   = 15,
        typEl       = 10000,
    )

    Tc = RealAxisSolver(inp)
    @test Tc == [8, 9]


    ### 9.TEST: Nb vDOS mu* - real axis, mu-update ###
    inp = arguments(
        a2f_file    = joinpath(@__DIR__, "Nb/Nb.a2f"),
        dos_file    = joinpath(@__DIR__, "Nb/Nb.dos"),
        Weep_file   = joinpath(@__DIR__, "Nb/Weep.dat"),
        flag_figure = 0,
        returnTc    = true,
        testMode    = true,
        outdir      = "./test/Nb/output/",
        cDOS_flag   = 0,
        mu_flag     = 1,
        ind_smear   = 15,
        typEl       = 10000,
    )

    Tc = RealAxisSolver(inp)
    @test Tc == [7, 8]


    ### 10.TEST: Nb vDOS W - real axis, fixed mu ###
    inp = arguments(
        a2f_file        = joinpath(@__DIR__, "Nb/Nb.a2f"),
        dos_file        = joinpath(@__DIR__, "Nb/Nb.dos"),
        Weep_file       = joinpath(@__DIR__, "Nb/Weep.dat"),
        flag_figure     = 0,
        returnTc        = true,
        testMode        = true,
        outdir          = "./test/Nb/output/",
        cDOS_flag       = 0,
        include_Weep    = 1,
        mu_flag         = 0,
        ind_smear       = 15,
        typEl           = 10000,
    )

    Tc = RealAxisSolver(inp)
    @test Tc == [7, 8]


    ### 11.TEST: Nb vDOS W - real axis, mu-update ###
    inp = arguments(
        a2f_file        = joinpath(@__DIR__, "Nb/Nb.a2f"),
        dos_file        = joinpath(@__DIR__, "Nb/Nb.dos"),
        Weep_file       = joinpath(@__DIR__, "Nb/Weep.dat"),
        flag_figure     = 0,
        returnTc        = true,
        testMode        = true,
        outdir          = "./test/Nb/output/",
        cDOS_flag       = 0,
        include_Weep    = 1,
        mu_flag         = 1,
        ind_smear       = 15,
        typEl           = 10000,
    )

    Tc = RealAxisSolver(inp)
    @test Tc == [7, 8]

end
