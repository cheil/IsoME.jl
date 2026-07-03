using IsoME


inp = arguments(
    a2f_file    = joinpath(@__DIR__, "Nb.a2F"),
    outdir      = joinpath(@__DIR__, "output"),
    cDOS_flag = 1,
    temps = [6],
    reOmega_c = 2000.0,
    domega = 0.4,
    mixing_beta = 0.5,
    mu = 0.38,
    typEl = 10000.0,
    conv_thr = 1e-3,
    ind_smear  = 15,
)

RealAxisSolver(inp)


