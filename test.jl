using IsoME

inp = arguments(
                    temps               = [20], 
                    a2f_file            = "/temp/spathd/MasterThesis/isotropic-me/input/TiN/TiN.a2f_qe_k24", 
                    outdir              = "/temp/spathd/MasterThesis/IsoME_Tests/",
                    dos_file            = "/temp/spathd/MasterThesis/isotropic-me/input/TiN/TiN_dos.dat",
                    Wen_file            = "/temp/spathd/MasterThesis/isotropic-me/input/TiN/DOS.dat",
                    Weep_file           = "/temp/spathd/MasterThesis/isotropic-me/input/TiN/Weep.dat",
                    TcSearchMode_flag   = 0,
                    cDOS_flag           = 0,
                    include_Weep        = 0
                )

EliashbergSolver(inp)


# fast testing
inp = arguments(
    temps = collect(10:20),
    cDOS_flag = 1,
    include_Weep = 0,
    a2f_file = "test/Nb/a2f_k48.dat",
    dos_file = "test/Nb/DOS_scf.dat",
    Weep_file = "test/Nb/Weep.dat",
    Wen_file = "test/Nb/Wenergies.dat",
    TcSearchMode_flag = 0,
    outdir = "/temp/spathd/MasterThesis/IsoME_Tests/output/",
    inp.conv_thr = 1e-2
    )


