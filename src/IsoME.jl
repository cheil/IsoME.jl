#
#
#  _                 __  __   ______ 
# | |               |  \/  | |  ____|
# | |  ___    ___   | \  / | | |__   
# | | / __|  / _ \  | |\/| | |  __|  
# | | \__ \ | (_) | | |  | | | |____ 
# |_| |___/  \___/  |_|  |_| |______|
#                                    
#                                    
#
#
# routine to solve the full-bandwidth isotropic Migdal-Eliashberg equations
# inspired by the EPW implementation
# 2023-10-17 - Christoph Heil


module IsoME


export EliashbergSolver, arguments, RealAxisSolver


using DelimitedFiles        # read in files   
using Interpolations        
using Plots, LaTeXStrings
using Trapz                 # integration
using LinearAlgebra         # 
using CSV                   # write to csv
using Printf                # Format console output, removable
using SparseIR              # Intermediate basis 
using LsqFit                # Curve fit in Tc search mode, removable
using Roots                 # Root finding
using Term                  # styled text terminal, removable
using Logging, LoggingExtras


### defining constants ###
const Ry2meV = 13605.662285137
const THz2meV = 4.13566553853599;
const kb = 0.08617333262; # meV/K

### Define input struct ###
# inputs Eliashberg Solver
@kwdef mutable struct arguments
    # Parameters
    temps::Vector{Float64}                  = [-1.0]        # concrete: Vector{Float64} is type-stable; Int entries auto-convert
    muc_AD::Float64                         = NaN           # NaN == "not assigned" (isnan check); NaN never a valid μ*
    imOmega_c::Float64                      = 7000.0
    muc_ME::Float64                         = NaN
    mu::Float64                             = NaN
    ef::Float64                             = NaN           # NaN == "auto-extract from DOS file header"
    efW::Float64                            = NaN           # NaN == "auto-extract from Weep file header"
    mixing_beta::Float64                    = NaN
    nItFullCoul::Int64                      = 10
    conv_thr::Float64                       = 1e-4
    minGap::Float64                         = 0.1
    N_it::Int64                             = 5000
    min_it::Int64                           = 10            # min iterations in eliashberg solver   
    encut::Float64                          = 5000.0      # symmetric outer energy cutoff: window [-encut, encut]
    shiftcut::Float64                       = 2000.0      # cutoff shift & Ne
    sparseSamplingTemp::Float64             = 2.0
    typEl::Float64                          = NaN
    flag_acon::Bool                         = false

    # interpolation
    itpStepSize::Vector{Int64}  = [1, 5, 50]
    itpBounds::Vector{Float64}  = [100.0, 500.0]

    # mode
    cDOS_flag::Int64    = 1
    include_Weep::Int64 = 0
    mu_flag::Int64      = 1
    broyden_flag::Int64 = 0     # self-consistency mixing: 0 = linear, 1 = Broyden (2nd method)
    broyden_mem::Int64  = 4     # Broyden history depth (number of stored iterations)

    # a2f input file
    a2f_file::String
    ind_smear::Int64    = -1            # -1 == auto (all header/footer/smear indices: -1 == auto-detect)
    nsmear::Int64       = -1
    nheader_a2f::Int64  = -1
    nfooter_a2f::Int64  = -1
    a2f_unit::String    = ""
    a2f_Nef::Float64    = NaN          # NaN == off; N(ε_F) used when α²F was computed. If set, α²F is
                                       # rescaled by dosef/a2f_Nef (dosef read from the DOS file)

    # dos input file
    dos_file::String    = ""
    nheader_dos::Int64  = -1
    nfooter_dos::Int64  = -1
    dos_unit::String    = ""
    spinDos::Int64      = 2

    # Weep input file
    Weep_file::String   = ""
    nheader_Weep::Int64 = -1
    nfooter_Weep::Int64 = -1
    Weep_unit::String   = ""
    Weep_col::Int64     = 3
    Wen_col::Int64      = 1
    Wen_file::String    = ""
    nheader_Wen::Int64  = -1
    nfooter_Wen::Int64  = -1
    Wen_unit::String    = ""

    # Output
    outdir::String              = pwd() 
    flag_figure::Int64          = 1
    flag_writeSelfEnergy::Int64 = 0
    material::String            = "Material"
    returnTc::Bool              = false
    testMode::Bool              = false

    # ----------- real axis inputs ----------- #
    # linear ω-grid
    reOmega_c::Float64          = 4000.0    # ω-grid cutoff
    domega::Float64             = 1.0       # ω-grid step size / meV
    dKernel::Float64            = 1.0       # (ω,ω') in Kernels
    # χ linear ω-grid
    reOmega_c_shift::Float64    = 25000.0   # ω-cutoff χ(ω)
    # ω'-chebyshev
    n_cheb::Int64               = 5000      # number of chebyshev points around poles in ω'-integration
    # ε-stepsize
    depsilon:: Int64            = 10


end


### include files ###
include("TcSearch.jl")
include("ReadIn.jl")
include("Interpolation.jl")
include("Mixing.jl")
include("AllenDynes.jl")
include("MuUpdate.jl")
include("WriteOutput.jl")
include("EliashbergEq.jl")
include("LinearKernels.jl")
include("wprimeGrid.jl")
include("realAxisSolver.jl")
include("realAxisEliashbergEq.jl")
include("Kernels.jl")
include("Acon.jl")


end