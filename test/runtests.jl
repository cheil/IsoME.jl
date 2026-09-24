# !!! Automatic test !!!
# Add further tests for encut, interpoaltion, ....

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
                    outdir      = joinpath(@__DIR__, "Nb", "output"),
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
                    outdir      = joinpath(@__DIR__, "Nb", "output"),
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
                    outdir      = joinpath(@__DIR__, "Nb", "output"),
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
    #     outdir      = joinpath(@__DIR__, "Nb", "output"),
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
        outdir      = joinpath(@__DIR__, "Nb", "output"),
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
        outdir      = joinpath(@__DIR__, "Nb", "output"),
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
        outdir      = joinpath(@__DIR__, "Nb", "output"),
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
        outdir      = joinpath(@__DIR__, "Nb", "output"),
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
        outdir      = joinpath(@__DIR__, "Nb", "output"),
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
        outdir          = joinpath(@__DIR__, "Nb", "output"),
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
        outdir          = joinpath(@__DIR__, "Nb", "output"),
        cDOS_flag       = 0,
        include_Weep    = 1,
        mu_flag         = 1,
        ind_smear       = 15,
        typEl           = 10000,
    )

    Tc = RealAxisSolver(inp)
    @test Tc == [7, 8]

    ############################################################
    # ------------------- sparse sampling -------------------- #
    ############################################################

    ### 12.TEST: Nb cDOS mu* at a single temperature, sparse Matsubara sampling ###
    # Sparse sampling kicks in below `sparseSamplingTemp` (default 2 K), which the Tc
    # searches above never reach - they all stop between 7 K and 10 K - so that branch of
    # eliashberg_eqn is otherwise never executed. Raising the threshold above the 5 K
    # solved here switches it on without paying for a sub-2 K run.

    # index set handed to the solver: sorted, starting at the first Matsubara frequency and
    # ending at the last one, and far smaller than the dense grid
    beta = 1 / (IsoME.kb * 5.0)
    omega_c = 7000.0
    M = ceil(Int, (omega_c / (pi * IsoME.kb * 5.0) - 1) / 2)
    ind_mat_freq = IsoME.initSparseSampling(beta, omega_c, M)
    @test issorted(ind_mat_freq)
    @test ind_mat_freq[1] == 1
    @test ind_mat_freq[end] == M + 1
    @test length(ind_mat_freq) < M + 1

    inp = arguments(
        a2f_file           = joinpath(@__DIR__, "Nb/Nb.a2f"),
        dos_file           = joinpath(@__DIR__, "Nb/Nb.dos"),
        Weep_file          = joinpath(@__DIR__, "Nb/Weep.dat"),
        flag_figure        = 0,
        returnTc           = true,
        testMode           = true,
        outdir             = joinpath(@__DIR__, "Nb", "output"),
        cDOS_flag          = 1,
        ind_smear          = 15,
        typEl              = 10000,
        temps              = [5.0],
        omega_c            = omega_c,
        sparseSamplingTemp = 10.0,
    )

    # 5 K is below Tc, so the gap survives: lower bound 5 K, no upper bound.
    # Same answer as the dense grid, which is what the sparse sampling has to reproduce.
    Tc = EliashbergSolver(inp)
    @test Tc[1] == 5.0
    @test isnan(Tc[2])

end
