# IsoME
IsoME is a state-of-the-art method to solve the isotropic Eliashberg equations. 
In addition to being the most accurate Eliashberg solver currently available, it also supports high-throughput 
calculations through a fast implementation based on the constant density of states (DOS) and Morel-Anderson pseudopotential approximations.

In its simplest form, IsoME requires only the Eliashberg spectral function ``\alpha^2 F`` as input to compute the superconducting critical temperature (Tc​). 
For more advanced calculations, input files containing the DOS, the screened Coulomb interaction (W), and the energy grid points of W are necessary.

The package is capable of solving the isotropic Eliashberg equations within any of the following approximations:
- constant Dos with Morel-Anderson pseudopotential ``\mu^*``
- variable Dos with Morel-Anderson pseudopotential ``\mu^*``
- constant Dos with screened Coulomb interaction ``W(\varepsilon,\varepsilon')``
- variable Dos with screened Coulomb interaction ``W(\varepsilon,\varepsilon')``

Furthermore, results based on the Allen-Dynes and Allen-Dynes-McMillan formulas are provided.  
There is also an option to save the self-energy components ``\Delta, Z, \chi`` at each temperature.  
An analytic continuation to the real frequency axis will be implemented in a future release.

## Installation
In order to run the code you need to [install](https://julialang.org/downloads/) julia 1.10 or higher.
IsoME.jl is not yet a registered package but it can be installed using the link to the github repository.
```julia-repl 
julia> using Pkg
julia> Pkg.add(path="https://github.com/cheil/IsoME.jl")
```

 ## Usage
After adding the package to your environment it can be loaded via
```julia-repl 
julia> using IsoME
```
To search for the Tc within the constant dos approximation using the Anderson-Pseuodopotential, only the path to the ``\alpha^2F``-file has to be provided. We recommend to set also the path to the output directory, beacuse per default the results will be written into the current working directory. The inputs are handed over by creating an instance of the custom data type `arguments()`.
```julia-repl 
julia> inp = arguments(
                a2f_file="Path to a2f-file", 
                outdir="Path to output directory"
                )
```
All other input values are optional and are contained within inp. 
The input arguments can be viewed by using a dot after inp, e.g.
```julia-repl 
julia> inp.omega_c
10000.0
```
gives the cutoff of the Matsubara summation in mEv. 
Note that some of the values in inp will only be set during execution. This is usually indicated by a -1 or "" for numeric and string values, respectively. 
E.g. if the ``\alpha^2F`` contains several columns with different smearing and the smearing index (ind_smear) is set to -1, the middle column will be used per default. 
The code should automatically recognize the format of the input files and remove the header and footer lines. If this fails, there is the option to specify these parameters manually. 
Furthermore, the unit and in the case of a vDOS or W calculation the fermi energy will be extracted from the header. For a detailed description have a look at the [Input Parameters](@ref) page.
Finally, the calculation can be started via
```julia-repl 
julia> EliashbergSolver(inp)
```
This creates a log-file and a file summarizing the results, as well as figures of the superconducting gap and the a2f-values. Saving the self-energy components has to be enabled explicitly.
Calculations within one of the other approximations can be started completely analog by specifying the paths to the dos/W-file and setting the corresponding flags.
```julia-repl   
julia> inp = arguments(
                a2f_file            = "Path to a2f-file", 
                outdir              = "Path to output directory",
                dos_file            = "Path to dos-file",
                Weep_file           = "Path to W-file",
                Wen_file            = "Path to W energies-file",
                cDOS_flag           = 0,
                include_Weep        = 1
                )
```
The Wen-file is optional if the Weep-file contains already the energy grid.

## Minimal example
A minimal example can be found at [Example files](https://github.com/cheil/IsoME.jl/tree/main/test/Nb).
If you have installed the package you should be able to run the file examples.jl in any of the four cases.
```console
~ $ julia --project= examples.jl
```


 ## Outline 
 ```@index
```

```@contents
```

