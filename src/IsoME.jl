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
# adapted from the EPW implementation
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
using Revise    # REMOVE before publishing


### defining constants ###
const Ry2meV = 13605.662285137
const THz2meV = 4.13566553853599;
const kb = 0.08617333262; # meV/K

### Define input struct ###
@kwdef mutable struct arguments
    # Parameters
    temps::Vector{Number} = [-1]        # change to type number?
    muc_AD::Float64 = 0.14
    omega_c::Float64 = 10000.0
    muc_ME::Float64 = -1
    mu::Float64 = -1  
    ef::Float64 = -1
    mixing_beta::Number = -1
    nItFullCoul::Number = 10
    conv_thr::Float64 = 1e-4
    N_it::Int64 = 5000             

    # mode
    cDOS_flag::Int64 = 1
    include_Weep::Int64 = 0
    TcSearchMode_flag::Int64 = 1
    mu_flag::Int64 = 1

    # a2f input file
    a2f_file::String
    ind_smear::Int64 = -1
    nsmear::Number   = -1
    nheader_a2f::Number = -1
    nfooter_a2f::Number = -1
    a2f_unit::String    = ""

    # dos input file
    dos_file::String = ""
    nheader_dos::Number = -1
    nfooter_dos::Number = -1
    dos_unit::String   = ""
    spinDos::Int64   = 2
    colFermi_dos::Int64 = 1

    # Weep input file, auto extraction if -1
    Weep_file::String = ""
    nheader_Weep::Number = -1
    nfooter_Weep::Number = -1
    Weep_unit::String = ""
    Weep_col::Number = 3
    Wen_file::String = ""
    nheader_Wen::Number = -1 
    nfooter_Wen::Number = -1
    Wen_unit::String   = ""

    # Output
    outdir::String = pwd()*"/output/" 
    flag_log::Int64 = 1
    flag_figure::Int64 = 1
    flag_outfile::Int64 = 1
    flag_writeSelfEnergy::Int64 = 0
    log_file::Any = ""
    material::String = "Material"

    # Restrict Weep
    Nrestrict::Number = -1
    wndRestrict::Vector{Number} = [-1]
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


end