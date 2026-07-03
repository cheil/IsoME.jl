using IsoME

inp = arguments(
    a2f_file    = joinpath(@__DIR__, "Al.alpha2F.dat"),
    outdir      = joinpath(@__DIR__, "output"),
    cDOS_flag = 1,
    temps = 0.1:0.1:2,
    real_c = 200,
    domega = 0.4,
    mixing_beta = 0.5,
    mu = 0.26,
    # muc_AD = 0.1,
    typEl = 10000,
    conv_thr = 1e-5,
    minGap = 0.01,
    ind_smear  = 1,
)

RealAxisSolver(inp)


