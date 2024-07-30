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
# Annotate type of global variables, e.g. x::Float64 = 1.0
# Remove revise before making module public
# print elapsed time 
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

### defining convergence parameters for Eliashberg solver ###
const N_it = 5000;        # max. number of iterations for the Eliashberg solver
const conv_thr = 1e-4;    # convergence threshold for Delta0 for solving the Eliashberg equations


### Define input struct ###
# Group into sub structs?
#= 
i.e. @kwdef struct a2f
    file::String = ""
    nheader::Number = nothing
end

@kwdef struct arguments 
    a2f::a2f = a2f(nheader = nothing)
end 
=#
@kwdef mutable struct arguments
    # Parameters
    temps::Vector{Number} = [-1]        # change to type number?
    muc::Float64 = 0.14
    omega_c::Float64 = 15000.0
    muc_ME::Float64 = muc / (1 + muc*log(200/omega_c))  
    mixing_beta::Number = -1
    nItFullCoul = 5             # if greater than 20 adapt termination criterion for min iterations !!

    # mode
    cDOS_flag::Int64 = 1
    TcSearchMode_flag::Int64 = 0
    mu_flag::Int64 = 1

    # a2f input file
    a2f_file::String
    ind_smear::Int64 = 1
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
    nheader_dosW::Number = -1 
    nfooter_dosW::Number = -1
    dosW_unit::String   = ""
    spinDosW::Number     = 1
    colFermi_dosW::Int64 = 0
    dosW_file::String = ""

    # Weep input file, auto extraction if -1
    Weep_file::String = ""
    nheader_Weep::Number = -1
    nfooter_Weep::Number = -1
    Weep_unit::String = ""
    include_Weep::Int64 = 0

    # Output
    outdir::String = "./output/"    # pwd() ??
    flag_log::Int64 = 1
    flag_figure::Int64 = 1
    flag_outfile::Int64 = 1
    log_file::Any = ""
    material::String = "Material"

    # Restrict Weep, REMOVE BEFORE MERGE !!!
    Nrestrict::Number = -1
    wndRestrict::Vector{Number} = [-1]
end


### include files ###
# only include files which are called directly here?
include("TcSearch.jl")
include("ReadIn.jl")
include("Interpolation.jl")
include("Mixing.jl")
include("AllenDynes.jl")
include("MuUpdate.jl")
include("FormatConsoleOutput.jl")
include("EliashbergEq.jl")



#=


### Start program ###
#try
#    dt = @elapsed include("IsoME_FindTc.jl")

#    print("Total Runtime: ", dt, " seconds\n")
#    if flag_log == 1
#        print(log_file, "Total Runtime: ", dt, " seconds\n")
    
#        # close & save
#        close(log_file)
#    end

#catch ex 

#    # write error to log file
#    showerror(log_file, ex, catch_backtrace())   
#    # close & save
#    close(log_file)

    # write error to console
#    showerror(stdout, ex, catch_backtrace())


#end


=#


end