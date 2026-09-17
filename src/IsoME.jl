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


using DelimitedFiles            # read in files   
using Interpolations        
using Plots, LaTeXStrings
using Trapz                     # integration
using LinearAlgebra             # 
using Printf                    # Format console output, also imported by other packages
using SparseIR                  # Intermediate basis, removable if we use a different approach for sparse sampling (Romans?)
using Logging, LoggingExtras    # also imported by other packages
using TOML



### defining constants ###
const Ry2meV = 13605.662285137
const THz2meV = 4.13566553853599;
const kb = 0.08617333262; # meV/K

### input validation (defines @checked_kwdef, used by the struct below) ###
include("InputValidation.jl")

### Define input struct ###
# inputs Eliashberg Solver
"""
    arguments(; a2f_file, kwargs...)

Every input of [`EliashbergSolver`](@ref) and [`RealAxisSolver`](@ref), collected in one
mutable struct. Only `a2f_file` is mandatory; all other fields have a default or are
inferred during the run.

Fields that are meant to be inferred carry a sentinel and are overwritten once their value
is known: `NaN` for real-valued fields, `-1` for integer fields (line counts, column and
smearing indices) and `""` for strings. Because the solver fills these in, an `arguments`
instance is *not* reusable — create a fresh one for every call.

All energies are in meV. Input files in meV, eV, THz, Ry or Ha are converted on read-in,
either from the unit in the file header or from the `*_unit` fields.

The approximation follows from two flags: `cDOS_flag` (1: constant DOS, 0: variable DOS)
and `include_Weep` (0: Morel-Anderson μ*, 1: static `W(ε,ε′)`). `RealAxisSolver` supports
cDOS+μ, vDOS+μ and vDOS+W; cDOS+W is imaginary-axis only and not recommended.

The fields are grouped as follows; the [Input](https://cheil.github.io/IsoME.jl/stable/Input/)
page of the documentation lists each one with its meaning and its default.

  * **Physical parameters** — `temps`, `omega_c`, `encut`, `mu`, `muc_AD`, `muc_ME`,
    `typEl`, `ef`, `efW`
  * **Mode** — `cDOS_flag`, `include_Weep`, `mu_flag`
  * **Convergence** — `conv_thr`, `N_it`, `min_it`, `minGap`, `mixing_beta`,
    `nItFullCoul`, `sparseSamplingTemp`
  * **Input files** — `a2f_file`, `dos_file`, `Weep_file`, `Wen_file`, together with their
    `nheader_*` / `nfooter_*` / `*_unit` / `*_col` fields, `ind_smear`, `nsmear`,
    `spinDos` and `a2f_Nef`
  * **ε-grid** — `itpBounds`, `itpStepSize` (imaginary axis), `depsilon` (real axis)
  * **Real axis** — `domega`, `n_cheb`
  * **Output** — `outdir`, `material`, `flag_figure`, `flag_writeSelfEnergy`, `flag_acon`,
    `returnTc`, `testMode`

Output is written to `outdir`, which defaults to `IsoME/` inside the current working
directory. An existing directory is never written into: a run counter is appended instead,
so repeated runs land in `IsoME_1/`, `IsoME_2/` and so on.

A wrong keyword name or a value of the wrong type is reported by field name, with a
suggestion, before the run starts. `mu`, `muc_AD`, `muc_ME`, `mixing_beta`, `conv_thr`,
`minGap`, `N_it` and `min_it` are additionally rejected when negative; `NaN` is unaffected,
so the "infer this during the run" sentinel keeps working.

# Examples
```julia
# Tc search, constant DOS with the default μ*_AD = 0.12
inp = arguments(a2f_file = "Nb.a2f", outdir = "Nb_cDOS")
Tc  = EliashbergSolver(inp)

# variable DOS with the full static Coulomb interaction, at two temperatures
inp = arguments(
    a2f_file     = "Nb.a2f",
    dos_file     = "Nb.dos",
    Weep_file    = "Weep.dat",
    cDOS_flag    = 0,
    include_Weep = 1,
    temps        = [5.0, 10.0],
    outdir       = "Nb_vDOS_W",
)
RealAxisSolver(inp)
```

See also [`EliashbergSolver`](@ref), [`RealAxisSolver`](@ref).
"""
@checked_kwdef mutable struct arguments
    # Parameters
    temps::Vector{Float64}                  = [-1.0]        # concrete: Vector{Float64} is type-stable; Int entries auto-convert
    muc_AD::Float64                         = NaN           # NaN == "not assigned" (isnan check); NaN never a valid μ*
    omega_c::Float64                        = 7000.0        # frequency cutoff, shared by both solvers (Matsubara / real axis)
    muc_ME::Float64                         = NaN
    mu::Float64                             = NaN
    ef::Float64                             = NaN           # NaN == "auto-extract from DOS file header"
    efW::Float64                            = NaN           # NaN == "auto-extract from Weep file header"
    mixing_beta::Float64                    = NaN
    nItFullCoul::Int64                      = 10
    conv_thr::Float64                       = 1e-4
    minGap::Float64                         = 0.1
    N_it::Int64                             = 5000
    min_it::Int64                           = 10            # min iterations, both solvers
    encut::Float64                          = 2000.0        # symmetric outer energy cutoff: window [-encut, encut]
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

    # a2f input file
    a2f_file::String
    ind_smear::Int64    = -1            # -1 == auto (all header/footer/smear indices: -1 == auto-detect)
    nsmear::Int64       = -1
    nheader_a2f::Int64  = -1
    nfooter_a2f::Int64  = -1
    a2f_unit::String    = ""
    a2f_Nef::Float64    = NaN          # NaN == off; N(ε_F) used when α²F was computed. If set, α²F is
                                       # rescaled by a2f_Nef/dosef (dosef read from the DOS file)

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
    outdir::String              = joinpath(pwd(), "IsoME")   # never the bare cwd: a run writes into its own directory
    flag_figure::Int64          = 1
    flag_writeSelfEnergy::Int64 = 0
    material::String            = "Material"
    returnTc::Bool              = false
    testMode::Bool              = false

    # ----------- real axis inputs ----------- #
    # linear ω-grid (cutoff: omega_c, see above)
    domega::Float64             = 1.0       # ω-grid step size / meV; also the kernel table step
    # ω'-chebyshev
    n_cheb::Int64               = 1000      # number of chebyshev points around poles in ω'-integration
    # ε-stepsize
    depsilon:: Int64            = 10

end


### include files ###
include("CurveFit.jl")
include("TcSearch.jl")
include("ReadIn.jl")
include("Interpolation.jl")
include("Mixing.jl")
include("AllenDynes.jl")
include("MuUpdate.jl")
include("WriteOutput.jl")
include("EliashbergEq.jl")
include("Kernels.jl")
include("wprimeGrid.jl")
include("realAxisSolver.jl")
include("realAxisEliashbergEq.jl")
include("Acon.jl")


### plotting defaults ###
"""
    setPlotDefaults()

Style shared by every figure IsoME writes, plus `show = false`: the solver saves its
figures to `outdir` and must never open a plot window, which would block a batch run.

Applied once when the module is loaded rather than from each plotting routine, so that the
style is defined in a single place and a run does not change the plotting defaults of the
session while it goes.
"""
function setPlotDefaults()
    default(
        show        = false,
        fontfamily  = "Computer Modern",
        linewidth   = 2,
        framestyle  = :box,
        label       = nothing,
        grid        = false,
    )
    return nothing
end

function __init__()
    setPlotDefaults()
    return nothing
end


### precompile workload ###
# build a dummy input structure, so that the (single, non-specialising) code path of the
# checked keyword constructor and of setproperty! ends up in the precompile cache
let
    inp = arguments(a2f_file = "")
    inp.temps = [1.0]
end


end