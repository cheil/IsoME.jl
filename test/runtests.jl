using IsoME
using Test


@testset "IsoME.jl" begin
    ### Nb ###
    inp = arguments(
                    temps = [17, 18, 19],
                    cDOS_flag = 1,
                    include_Weep = 0,
                    a2f_file = "test/Nb/a2f_k48.dat",
                    dos_file = "test/Nb/DOS_scf.dat",
                    Weep_file = "test/Nb/Weep.dat",
                    Wen_file = "test/Nb/Wenergies.dat",
                    flag_log = 0,
                    flag_figure = 0,
                    flag_outfile = 0
                    )

    Tc = EliashbergSolver(inp, true)
    @test Tc == 18
       
    ### USE CORRECT Tc VALUES !!! ###
    inp.cDOS_flag = 0
    inp.temps = [23,24,25]
    Tc = EliashbergSolver(inp, true)
    @test Tc == 26

    inp.cDOS_flag = 1
    inp.include_Weep = 1
    inp.temps = [26,27,28]
    Tc = EliashbergSolver(inp, true)
    @test Tc == 27

    inp.cDOS_flag = 0
    inp.temps = [23,24,25]
    Tc = EliashbergSolver(inp, true)
    @test Tc == 24
    
end


