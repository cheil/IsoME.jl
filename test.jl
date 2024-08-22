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


