#
#
#  _                 __  __   ______ 
# (_)               |  \/  | |  ____|
#  _   ___    ___   | \  / | | |__   
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

# ToDo:
# Remove revise before making module public
# error handle


export EliashbergSolver, arguments


using DelimitedFiles        # read in files
using Interpolations    
using Plots, LaTeXStrings
using Trapz                 # integration
using LinearAlgebra         # 
using CSV, DataFrames       # write to csv
using Printf                # Format console output
using SparseIR              # Intermediate basis 
using LsqFit                # Fitting 
using Roots                 # Root finding
using Term                  # styled text terminal
using Nevanlinna            # analytic continuation
using Logging, LoggingExtras
using Revise                # !!! REMOVE before publishing


### defining constants ###
const Ry2meV = 13605.662285137
const THz2meV = 4.13566553853599;
const kb = 0.08617333262; # meV/K


### Define input struct ###
@kwdef mutable struct arguments
    # Parameters
    temps::Vector{Number}   = [-1]        
    muc_AD::Float64         = -1
    omega_c::Float64        = 7000.0
    muc_ME::Float64         = -1
    mu::Float64             = -1  
    ef::Float64             = -1
    efW::Float64            = -1
    mixing_beta::Number     = -1
    nItFullCoul::Number     = 10
    conv_thr::Float64       = 1e-4
    minGap::Float64         = 0.1
    N_it::Int64             = 5000   
    encut::Float64           = -1        # outer cutoff energies
    shiftcut::Float64          = 2000      # cutoff shift & Ne
    sparseSamplingTemp::Float64 = 5  
    typEl::Float64          = -1 
    
    # interpolation
    itpStepSize::Vector{Int64}  = [1, 5, 50]
    itpBounds::Vector{Float64}  = [50, 200]

    # mode
    cDOS_flag::Int64    = 1
    include_Weep::Int64 = 0
    mu_flag::Int64      = 1

    # a2f input file
    a2f_file::String
    ind_smear::Int64    = -1
    nsmear::Int64       = -1
    nheader_a2f::Int64  = -1
    nfooter_a2f::Int64  = -1
    a2f_unit::String    = ""

    # dos input file
    dos_file::String    = ""
    nheader_dos::Number = -1
    nfooter_dos::Number = -1
    dos_unit::String    = ""
    spinDos::Int64      = 2

    # Weep input file
    Weep_file::String   = ""
    nheader_Weep::Int64 = -1
    nfooter_Weep::Int64 = -1
    Weep_unit::String   = ""
    Weep_col::Int64     = 3
    Wen_file::String    = ""
    nheader_Wen::Int64  = -1 
    nfooter_Wen::Int64  = -1
    Wen_unit::String    = ""

    # Output
    outdir::String      = pwd() 
    flag_figure::Int64  = 1
    flag_writeSelfEnergy::Int64 = 0
    material::String    = "Material"
    returnTc::Bool     = false

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
include("Integration.jl")
include("Acon.jl")


end