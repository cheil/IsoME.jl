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


export my_f, a
#    findTc


# use packages
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


a = 1

my_f(x,y) = 2x+3y

# include files
# rename w/o IsoME?
# only include files which are called directly here?
#include("ReadIn.jl")
#include("Interpolation.jl")
#include("Mixing.jl")
#include("AllenDynes.jl")
#include("MuUpdate.jl")
#include("FormatConsoleOutput.jl")
#include("EliashbergEq.jl")
#include("FindTc.jl")


#=
### Defining some defaults ###s
# Parameters
(@isdefined muc) || (muc = 0.1)
(@isdefined temps) || (temps = collect(5.0:5.0:200.0))
(@isdefined omega_c) || (omega_c = 20.0 * a2f_omega_raw[end])

# mode
(@isdefined cDOS_flag) || (cDOS_flag = 0)
(@isdefined TcSearchMode_flag) || (TcSearchMode_flag = 0)
(@isdefined mu_flag) || (mu_flag = 1)

# a2f input file, auto extraction if nothing
(@isdefined ind_smear) || (ind_smear = 1)
(@isdefined nsmear) || (nsmear = nothing)
(@isdefined nheader_a2f) || (nheader_a2f = nothing)
(@isdefined nfooter_a2f) || (nfooter_a2f = nothing)
(@isdefined a2f_unit) || (a2f_unit = nothing)

# dos input file, auto extraction if nothing
(@isdefined nheader_dos) || (nheader_dos = nothing)
(@isdefined nfooter_dos) || (nfooter_dos = nothing)
(@isdefined dos_unit) || (dos_unit = nothing)
(@isdefined spinDos) || (spinDos = 2)
(@isdefined colFermi_dos) || (colFermi_dos = 1)
(@isdefined nheader_dosW) || (nheader_dosW = nothing)
(@isdefined nfooter_dosW) || (nfooter_dosW = nothing)
(@isdefined dosW_unit) || (dosW_unit = nothing)
(@isdefined spinDosW) || (spinDosW = 1)
(@isdefined colFermi_dosW) || (colFermi_dosW = 0)

# Weep input file, auto extraction if nothing
(@isdefined nheader_Weep) || (nheader_Weep = nothing)
(@isdefined nfooter_Weep) || (nfooter_Weep = nothing)
(@isdefined Weep_unit) || (Weep_unit = nothing)
(@isdefined include_Weep) || (include_Weep = 1)

# Output
(@isdefined outdir) || (outdir = "./output/" )
(@isdefined flag_log) || (flag_log = 1)
(@isdefined flag_figure) || (flag_figure = 1)
(@isdefined flag_outfile) || (flag_outfile = 1)

# Restrict Weep, REMOVE BEFORE MERGE !!!
(@isdefined Nrestrict) || (Nrestrict = nothing)
(@isdefined wndRestrict) || (wndRestrict = nothing)




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
