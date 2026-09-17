# IsoME
IsoME solves the isotropic Eliashberg equations on either the imaginary or the real frequency axis, through the functions [`EliashbergSolver()`](@ref "Matsubara Solver") and [`RealAxisSolver()`](@ref "Real Axis Solver").
Both solvers collect their input parameters in the custom `arguments()` type, and both are driven by the same input files.

In its simplest form, IsoME requires only the Eliashberg spectral function ``\alpha^2 F`` to estimate the superconducting critical temperature ``\mathrm{T}_C``.
For more advanced calculations, input files containing the electronic density of states (DOS) and the screened Coulomb interaction ``W`` can be supplied, which allows the Eliashberg equations to be solved in any of the following approximations:
- cDOS``+\mu``: constant DOS with Morel-Anderson pseudopotential ``\mu^*``
- vDOS``+\mu``: variable DOS with Morel-Anderson pseudopotential ``\mu^*``
- vDOS``+W``: variable DOS with screened Coulomb interaction ``W(\varepsilon,\varepsilon')``

A fourth combination, the static Coulomb interaction ``W(\varepsilon,\varepsilon')`` together with the constant-DOS approximation, is also implemented on the imaginary axis.
However, this variant is not recommended, as it requires the same input data as the vDOS``+W`` approach while being less rigorous and offering no notable computational advantage.

Alongside the Eliashberg results, estimates from the Allen-Dynes-McMillan formula and from its machine-learning-corrected variant are always reported, and the self-energy components ``\Delta, Z, \chi, \phi`` can be saved at each temperature.

## The two solvers
Imaginary-axis calculations are started with [`EliashbergSolver()`](@ref "Matsubara Solver"), which is based on the [IsoME](https://doi.org/10.1016/j.cpc.2025.109720) paper.
The imaginary-axis solver is fast and robust, which makes it the natural choice for ``\mathrm{T}_C`` searches; its constant-DOS mode is particularly suited to high-throughput calculations.
To obtain spectral quantities, the converged solution can be analytically continued to real frequencies with the built-in Padé continuation.

Direct real-axis calculations are started with [`RealAxisSolver()`](@ref "Real Axis Solver"), which follows the real-axis formulation of [Simon *et al.*](https://doi.org/10.48550/arXiv.2603.18199) and supports the same three approximations.
It solves for the complex gap, renormalization and shift functions directly on a real-frequency grid, and is therefore the tool of choice whenever real-frequency self-energy components are needed without relying on an analytic continuation.
The price is a numerically more demanding calculation, so the real-axis solver is more sensitive to its grids and cutoffs than the imaginary-axis solver — see the [FAQ](@ref) for common pitfalls.

## Installation
To run the code, you need to [install](https://julialang.org/downloads/) Julia 1.10 or higher.
IsoME.jl is a registered package and can be installed using the Julia package manager.
```julia-repl
julia> using Pkg
julia> Pkg.add("IsoME")
```

## Usage
After adding the package to your environment, it can be loaded via
```julia-repl
julia> using IsoME
```
To search for ``\mathrm{T}_C`` within the constant-DOS approximation using the Morel-Anderson pseudopotential, only the path to the ``\alpha^2F`` file has to be provided.
We also recommend setting the output directory explicitly; otherwise results are written to `IsoME/` inside the current working directory.
An existing directory is never written into — a counter is appended instead, so repeated runs land in `IsoME_1/`, `IsoME_2/`, and so on.
Inputs are collected by creating an instance of `arguments()`.
```julia-repl
julia> inp = arguments(
                a2f_file="Path to a2f-file",
                outdir="Path to output directory"
                )
```
All other input values are optional and are contained within `arguments()`.
The inputs can be viewed using dot notation, e.g.
```julia-repl
julia> inp.omega_c
7000.0
```
gives the frequency cutoff in meV — the Matsubara cutoff on the imaginary axis, the ``\omega``-grid cutoff on the real axis.
Values that are determined during execution are marked with `NaN` for real-valued fields, `-1` for integer fields and `""` for strings.
For example, if the ``\alpha^2F`` file contains several smearing columns and `ind_smear` is left unset, the middle column is used by default.
In the same spirit, IsoME detects the format of the input files automatically: it strips header and footer lines and reads the units from the header, together with the Fermi energy for vDOS or ``W`` calculations.
For a detailed description, see the [Input](@ref) page.
Finally, the calculation can be started via
```julia-repl
julia> EliashbergSolver(inp)
```
The solver writes a log file, a result summary, an overview of the inputs, and, if enabled, figures of the superconducting gap and the ``\alpha^2F`` values.
The self-energy components are saved and plotted (`.dat` + `.png`, in `outdir/SelfEnergy/`) only when `flag_writeSelfEnergy = 1`.

Other approximations are selected by providing the DOS or ``W`` file paths and setting the corresponding flags.
```julia-repl
julia> inp = arguments(
                a2f_file            = "Path to a2f-file",
                outdir              = "Path to output directory",
                dos_file            = "Path to dos-file",
                Weep_file           = "Path to W-file",
                cDOS_flag           = 0,
                include_Weep        = 1
                )
```
If the energy grid is not contained in the `Weep_file`, it can be provided through an additional file: `Wen_file = "Path to W energies-file"`.

The same input structure drives the real-axis solver, which is started with `RealAxisSolver()` instead:
```julia-repl
julia> inp = arguments(
                a2f_file            = "Path to a2f-file",
                outdir              = "Path to output directory"
                )

julia> RealAxisSolver(inp)
```
Here, too, the approximation follows from the input files and flags, and the ``\mu`` and ``\mu^*`` values are handled exactly as on the imaginary axis (cDOS``+W`` is available on the imaginary axis only).
The parameters that control its grids are described under [Real-axis inputs](@ref "Real-axis inputs"), and their effect is explained on the [Real Axis Solver](@ref) page.

## Minimal example
A minimal example can be found at [Example files](https://github.com/cheil/IsoME.jl/tree/main/test/Nb).
If you have installed the package, you should be able to run `examples.jl` in any of the supported approximations.
```console
~ $ julia examples.jl
```
If you have installed IsoME into a separate project environment, which is the recommended way, specify the path to the environment via
```console
~ $ julia --project=/path/to/environment/ examples.jl
```
For more information, please refer to [Best practices](@ref).

