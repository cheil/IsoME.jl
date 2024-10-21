using IsoME

#=
# Al
inp = arguments(
    temps               = [5], 
    material            = "Al",
    a2f_file            = "/temp/spathd/MasterThesis/isotropic-me/input/Al/Al.a2F", 
    outdir              = "/temp/spathd/MasterThesis/IsoME_Tests/Al"*"_Weep",
    dos_file            = "/temp/spathd/MasterThesis/isotropic-me/input/Al/Al.dos",
    Wen_file            = "/temp/spathd/MasterThesis/isotropic-me/input/Al/DOS.dat",
    Weep_file           = "/temp/spathd/MasterThesis/isotropic-me/input/Al/Weep.dat",
    TcSearchMode_flag   = 1,
    cDOS_flag           = 0,
    include_Weep        = 1,
    omega_c             = 10000,
    mu_flag             = 1
)

EliashbergSolver(inp)


# Al
inp = arguments(
    temps               = [5], 
    material            = "Al",
    a2f_file            = "/temp/spathd/MasterThesis/isotropic-me/input/Al/Al.a2F", 
    outdir              = "/temp/spathd/MasterThesis/IsoME_Tests/Al"*"_muc",
    dos_file            = "/temp/spathd/MasterThesis/isotropic-me/input/Al/Al.dos",
    Wen_file            = "/temp/spathd/MasterThesis/isotropic-me/input/Al/DOS.dat",
    Weep_file           = "/temp/spathd/MasterThesis/isotropic-me/input/Al/Weep.dat",
    TcSearchMode_flag   = 1,
    cDOS_flag           = 0,
    include_Weep        = 0,
    omega_c             = 10000,
    mu_flag             = 1
)

EliashbergSolver(inp)



# TiN
inp = arguments(
    temps               = [20], 
    a2f_file            = "/temp/spathd/MasterThesis/isotropic-me/input/TiN/TiN.a2f_qe_k24", 
    outdir              = "/temp/spathd/MasterThesis/tests/TiN_temp",
    dos_file            = "/temp/spathd/MasterThesis/isotropic-me/input/TiN/TiN_dos.dat",
    Wen_file            = "/temp/spathd/MasterThesis/isotropic-me/input/TiN/DOS.dat",
    Weep_file           = "/temp/spathd/MasterThesis/isotropic-me/input/TiN/Weep.dat",
    TcSearchMode_flag   = 0,
    cDOS_flag           = 0,
    include_Weep        = 0,
    omega_c             = 10000,
    mu_flag             = 1,
    fsthick =           10000
)
EliashbergSolver(inp)


inp = arguments(
    temps               = [20], 
    a2f_file            = "/temp/spathd/MasterThesis/isotropic-me/input/TiN/TiN.a2f_qe_k24", 
    outdir              = "/temp/spathd/MasterThesis/IsoME_Tests/TiN_muc",
    dos_file            = "/temp/spathd/MasterThesis/isotropic-me/input/TiN/TiN_dos.dat",
    Wen_file            = "/temp/spathd/MasterThesis/isotropic-me/input/TiN/DOS.dat",
    Weep_file           = "/temp/spathd/MasterThesis/isotropic-me/input/TiN/Weep.dat",
    TcSearchMode_flag   = 1,
    cDOS_flag           = 0,
    include_Weep        = 0,
    omega_c             = 10000,
    mu_flag             = 1
)
EliashbergSolver(inp)
=#

#=
# NbN
inp = arguments(
    material="NbN",
    temps=[20],
    a2f_file="/temp/spathd/MasterThesis/isotropic-me/input/NbN/NbN.a2f_epw_kq24_degauss02",
    outdir="/temp/spathd/MasterThesis/IsoME_Tests/NbN_" * string(indSmear),
    dos_file="/temp/spathd/MasterThesis/isotropic-me/input/NbN/NbN.dos_qe.dat",
    Wen_file="/temp/spathd/MasterThesis/isotropic-me/input/NbN/DOS.dat",
    Weep_file="/temp/spathd/MasterThesis/isotropic-me/input/NbN/Weep.dat",
    TcSearchMode_flag=1,
    cDOS_flag=0,
    include_Weep=1,
    omega_c=10000,
    mu_flag=1,
    ind_smear=indSmear,
)
EliashbergSolver(inp)
=#
#=
inp = arguments(
    temps               = [20], 
    a2f_file            = "/temp/spathd/MasterThesis/isotropic-me/input/NbN/NbN.a2f_epw_kq24_degauss02", 
    outdir              = "/temp/spathd/MasterThesis/IsoME_Tests/NbN_muc",
    dos_file            = "/temp/spathd/MasterThesis/isotropic-me/input/NbN/NbN.dos_qe.dat",
    Wen_file            = "/temp/spathd/MasterThesis/isotropic-me/input/NbN/DOS.dat",
    Weep_file           = "/temp/spathd/MasterThesis/isotropic-me/input/NbN/Weep.dat",
    TcSearchMode_flag   = 1,
    cDOS_flag           = 0,
    include_Weep        = 0,
    omega_c             = 10000,
    mu_flag             = 1
)
EliashbergSolver(inp)
=#


#=
# Nb 
inp = arguments(
    temps = [10],
    cDOS_flag = 0,
    include_Weep = 0,
    a2f_file = "test/Nb/a2f_k48.dat",
    dos_file = "test/Nb/DOS_scf.dat",
    Weep_file = "test/Nb/Weep.dat",
    Wen_file = "test/Nb/Wenergies.dat",
    TcSearchMode_flag = 0,
    outdir = "/temp/spathd/MasterThesis/tests/temp/",
    mu = 0.1735,
    )


# fast testing
inp = arguments(
    temps = [10],
    cDOS_flag = 1,
    include_Weep = 0,
    a2f_file = "test/Nb/a2f_k48.dat",
    dos_file = "test/Nb/DOS_scf.dat",
    Weep_file = "test/Nb/Weep.dat",
    Wen_file = "test/Nb/Wenergies.dat",
    TcSearchMode_flag = 0,
    outdir = "/temp/spathd/MasterThesis/tests/temp/",
    conv_thr = 1e-2
    )
=#


# Nb_Antonio
#=
inp = arguments(
    temps = [12,13,14,15,16],
    cDOS_flag = 1,
    include_Weep = 0,
    a2f_file = "/temp/spathd/materials_IsoME_paper/Nb_Antonio/Nb_Antonio.a2F",
    dos_file = "/temp/spathd/materials_IsoME_paper/Nb_Antonio/Nb_Antonio.dos",
    Weep_file = "/temp/spathd/materials_IsoME_paper/Nb_Antonio/Weep.dat",
    Wen_file = "/temp/spathd/materials_IsoME_paper/Nb_Antonio/Wen.dat",
    TcSearchMode_flag = 0,
    outdir = "/temp/spathd/MasterThesis/tests/temp/",
    mu = -1,
    )

inp = arguments(
    temps = [12,13,14,15,16],
    cDOS_flag = 1,
    include_Weep = 0,
    a2f_file = "/temp/spathd/MasterThesis/materials/Nb_Antonio/elias_ph.in",
    dos_file = "/temp/spathd/MasterThesis/materials/Nb_Antonio/dos.dat",
    Weep_file = "/temp/spathd/MasterThesis/materials/Nb_Antonio/KC.OUT",
    Wen_file = "/temp/spathd/MasterThesis/materials/Nb_Antonio/Wen.dat",
    TcSearchMode_flag = 0,
    outdir = "/temp/spathd/MasterThesis/tests/temp/",
    mu = -1,
    )
=#