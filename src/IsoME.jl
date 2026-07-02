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
@kwdef mutable struct arguments{T<:Number}
    # Parameters
    temps::Union{Vector{T}, Nothing}   = nothing
    muc_AD::Union{Float64, Nothing}         = nothing
    imOmega_c::Float64      = 7000.0
    muc_ME::Union{Float64, Nothing}         = nothing
    mu::Union{Float64, Nothing}             = nothing
    ef::Union{Float64, Nothing}  = nothing
    efW::Union{Float64, Nothing} = nothing
    mixing_beta::Union{Float64, Nothing}     = nothing
    nItFullCoul::Int64     = 10
    conv_thr::Float64       = 1e-4
    minGap::Float64         = 0.1
    N_it::Int64             = 5000
    min_it::Int64           = 10            # min iterations in eliashberg solver   
    encut::Union{Float64, Vector{Float64}, Nothing} = nothing      # outer cutoff energies
    shiftcut::Float64          = 2000.0      # cutoff shift & Ne
    sparseSamplingTemp::Float64 = 2.0
    typEl::Union{Float64, Nothing} = nothing
    flag_acon::Bool         = false   
    plot_flag::Bool         = false  
    
    # interpolation
    itpStepSize::Vector{Int64}  = [1, 5, 50]
    itpBounds::Vector{Float64}  = [100, 500]

    # mode
    cDOS_flag::Int64    = 1
    include_Weep::Int64 = 0
    mu_flag::Int64      = 1

    # a2f input file
    a2f_file::String
    ind_smear::Union{Int64, Nothing}    = nothing
    nsmear::Union{Int64, Nothing}       = nothing
    nheader_a2f::Union{Int64, Nothing}  = nothing
    nfooter_a2f::Union{Int64, Nothing}  = nothing
    a2f_unit::String    = ""

    # dos input file
    dos_file::String    = ""
    nheader_dos::Union{Int64, Nothing}  = nothing
    nfooter_dos::Union{Int64, Nothing}  = nothing
    dos_unit::String    = ""
    spinDos::Int64      = 2

    # Weep input file
    Weep_file::String   = ""
    nheader_Weep::Union{Int64, Nothing} = nothing
    nfooter_Weep::Union{Int64, Nothing} = nothing
    Weep_unit::String   = ""
    Weep_col::Int64     = 3
    Wen_col::Int64      = 1
    Wen_file::String    = ""
    nheader_Wen::Union{Int64, Nothing}  = nothing
    nfooter_Wen::Union{Int64, Nothing}  = nothing
    Wen_unit::String    = ""

    # Output
    outdir::String      = pwd() 
    flag_figure::Int64  = 1
    flag_writeSelfEnergy::Int64 = 0
    material::String    = "Material"
    returnTc::Bool      = false
    testMode::Bool      = false

    # ----------- real axis inputs ----------- #
    # ω-grid
    reOmega_c::Float64          = 2000.0  # ω-grid cutoff
    numReal_c::Int64            = 4000  # num of ω-points
    # χ ω-grid
    reOmega_c_shift::Float64    = 15000.0 # ω-cutoff χ(ω)
    numReal_c_shift::Int64      = 10000 # num of ω-points χ
    # Ω-chebyshev
    n_cheb::Int64               = 5000  # number of chebyshev points in Ω-integration
    # ω'-grid
    num_wp1::Int64              = 1000  # inner chebyshev grid
    num_wp2::Int64              = 5000  # outer chebyshev grid
    wp_max::Float64             = 2.0     # Location outer chebyshev grid
    # ε-stepsize
    depsilon:: Int64            = 10


end

# Default the (temps-only) type parameter to Float64 when it can't be inferred
# from the keyword arguments (e.g. when temps is left as nothing). Users can
# still force integer temperatures via arguments{Int}(; ...).
arguments(; kwargs...) = arguments{Float64}(; kwargs...)


### include files ###
include("TcSearch.jl")
include("ReadIn.jl")
include("Interpolation.jl")
include("Mixing.jl")
include("AllenDynes.jl")
include("MuUpdate.jl")
include("WriteOutput.jl")
include("EliashbergEq.jl")
include("realAxisSolver.jl")
include("realAxisEliashbergEq.jl")
include("Kernels.jl")
include("Acon.jl")


end