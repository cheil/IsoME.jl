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
# Transfer ReadMe to github
# Annotate type of global variables, e.g. x::Float64 = 1.0
# Default values struct
# Remove revise before making module public

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


### defining constants ###
const Ry2meV = 13605.662285137
const THz2meV = 4.13566553853599;
const kb = 0.08617333262; # meV/K

### defining convergence parameters for Eliashberg solver ###
const N_it = 5000;        # max. number of iterations for the Eliashberg solver
const conv_thr = 1e-4;    # convergence threshold for Delta0 for solving the Eliashberg equations


### Define input struct ###
@kwdef struct arguments
    # Parameters
    temps::Vector{Float64}
    muc::Float64 = 0.16
    omega_c::Float64 = 15000

    # mode
    cDOS_flag::Int64 = 1
    TcSearchMode_flag::Int64 = 0
    mu_flag::Int64 = 1

    # a2f input file
    a2f_file::String
    ind_smear::Int64 = 1
    nsmear::Int64   = Int64
    nheader_a2f::Int64 = Int64
    nfooter_a2f::Int64 = Int64
    a2f_unit::Int64    = Int64

    # dos input file
    dos_file::String = ""
    nheader_dos::Int64 = Int64 
    nfooter_dos::Int64 = Int64 
    dos_unit::Int64   = Int64
    spinDos::Int64   = 2
    colFermi_dos::Int64 = 1
    nheader_dosW::Int64 = Int64 
    nfooter_dosW::Int64 = Int64
    dosW_unit::Int64   = Int64
    spinDosW::Int64     = 1
    colFermi_dosW::Int64 = 0
    dosW_file::String = ""

    # Weep input file, auto extraction if nothing
    Weep_file::String = ""
    nheader_Weep::Int64 = Int64
    nfooter_Weep::Int64 = Int64
    Weep_unit::Int64 = Int64
    include_Weep::Int64 = 0

    # Output
    outdir::String = "./output/"
    flag_log::Int64 = 1
    flag_figure::Int64 = 1
    flag_outfile::Int64 = 1 

    # Restrict Weep, REMOVE BEFORE MERGE !!!
    Nrestrict::Int64 = Int64
    wndRestrict::Vector{Float64} = Vector{Int64}()
end

### include files ###
# only include files which are called directly here?
include("ReadIn.jl")
include("Interpolation.jl")
include("Mixing.jl")
include("AllenDynes.jl")
include("MuUpdate.jl")
include("FormatConsoleOutput.jl")
include("EliashbergEq.jl")
include("TcSearch.jl")


#=
### open log_file ###
#if flag_log == 1
#    log_file = open(outdir * "/log.txt", "w")
#end

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