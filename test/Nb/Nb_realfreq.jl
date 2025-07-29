using IsoME


inp = arguments(
    a2f_file    = joinpath(@__DIR__, "Nb.a2F"),
    outdir      = joinpath(@__DIR__, "output"),
    cDOS_flag = 1,
    temps = 1:1:10,
    real_c = 200,
    numReal_c = 5000,
    mixing_beta = 0.5,
    mu = 0.38,
    typEl = 10000,
    conv_thr = 1e-5,
    minGap = 0.01,
    ind_smear  = 15,
)

RealAxisSolver(inp)


