# IsoME
IsoME solves the isotropic Eliashberg equations on either the imaginary or the real frequency axis.
Imaginary-axis calculations are started with `EliashbergSolver()`, while direct real-axis calculations are started with `RealAxisSolver()`.
Both solvers use the custom `arguments()` type to collect input parameters.

The imaginary-axis solver `EliashbergSolver()` is based on the [IsoME](https://doi.org/10.1016/j.cpc.2025.109720) paper.
It provides an accurate implementation of the isotropic Migdal-Eliashberg equations and is also suited for high-throughput calculations through a fast constant-density-of-states (cDOS) mode with a Morel-Anderson pseudopotential.

In its simplest form, IsoME requires only the Eliashberg spectral function ``\alpha^2 F`` to estimate the superconducting critical temperature ``T_c``.
For more advanced calculations, input files containing the electronic density of states (DOS) and the screened Coulomb interaction ``W`` can be supplied.

The package is capable of solving the isotropic Eliashberg equations within any of the following approximations:
- cDOS``+\mu``: constant DOS with Morel-Anderson pseudopotential ``\mu^*``
- vDOS``+\mu``: variable DOS with Morel-Anderson pseudopotential ``\mu^*``
- vDOS``+W``: variable DOS with screened Coulomb interaction ``W(\varepsilon,\varepsilon')``

Furthermore, results based on the Allen-Dynes and Allen-Dynes-McMillan formulas are provided.
In principle, a fourth level of approximation that combines the static Coulomb interaction ``W(\varepsilon,\varepsilon')`` with the constant DOS approximation exists.
However, this variant is not recommended, as it requires the same input data as the vDOS+W approach while being less rigorous and offering no notable computational advantage.
There is also an option to save the self-energy components ``\Delta, Z, \chi, \phi`` at each temperature.

The imaginary-axis solution can also be analytically continued to the real frequency axis with the implemented Pade approximation.
This is often the most convenient way to inspect spectral quantities after an imaginary-axis run.

The direct real-axis solver `RealAxisSolver()` follows the real-axis formulation of [Simon](https://doi.org/10.48550/arXiv.2603.18199).
It solves for the frequency-dependent superconducting gap and renormalization functions directly on a real-frequency grid.
The real-axis solver currently supports the ``\mu`` approximation for cDOS and vDOS calculations; screened-Coulomb ``W`` will be implemented in the future.
It is in particular useful when real-frequency self-energy components are desired without relying on analytic continuation.

## Installation
In order to run the code you need to [install](https://julialang.org/downloads/) Julia 1.10 or higher.
IsoME.jl is a registered package and it can be installed using the Julia package manager.
```julia-repl
julia> using Pkg
julia> Pkg.add("IsoME")
```

## Usage
After adding the package to your environment it can be loaded via
```julia-repl
julia> using IsoME
```
To search for ``T_c`` within the constant-DOS approximation using the Morel-Anderson pseudopotential, only the path to the ``\alpha^2F`` file has to be provided.
We also recommend setting the output directory explicitly; otherwise, results are written to the current working directory.
Inputs are collected by creating an instance of `arguments()`.
```julia-repl
julia> inp = arguments(
                a2f_file="Path to a2f-file",
                outdir="Path to output directory"
                )
```
All other input values are optional and are contained within `arguments()`.
The inputs can be viewed using the dot notation, e.g.
```julia-repl
julia> inp.imOmega_c
7000.0
```
gives the cutoff of the Matsubara summation in meV.
Some values in `inp` are determined during execution; these are usually marked by `-1` for numeric values or `""` for strings.
For example, if the ``\alpha^2F`` file contains several smearing columns and `ind_smear` is left at `-1`, the middle column is used by default.
IsoME attempts to recognize input-file formatting automatically and removes header and footer lines.
It also extracts units and, for vDOS or ``W`` calculations, the Fermi energy from the file headers when possible.
For a detailed description, see the [Input](@ref) page.
Finally, the calculation can be started via
```julia-repl
julia> EliashbergSolver(inp)
```
The solver writes a log file, a result summary, an overview of the inputs, and, if enabled, figures of the superconducting gap and the ``\alpha^2F`` values.
Self-energy components are saved only when `flag_writeSelfEnergy = 1`.
Other approximations can be selected by providing the DOS or ``W`` file paths and setting the corresponding flags.
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

For solving the Eliashberg equations directly on the real axis, use `RealAxisSolver()` instead:
```julia-repl
julia> inp = arguments(
                a2f_file            = "Path to a2f-file",
                outdir              = "Path to output directory",
                )

julia> RealAxisSolver(inp)
```

## Minimal example
A minimal example can be found at [Example files](https://github.com/cheil/IsoME.jl/tree/main/test/Nb).
If you have installed the package, you should be able to run `examples.jl` in any of the supported approximations.
```console
~ $ julia examples.jl
```
If you have installed IsoME into a separate project environment, specify the path to the environment via
```console
~ $ julia --project=/path/to/environment/ examples.jl
```
For more information please refer to [Best practices](@ref).

