# IsoME
IsoME.jl offers a quick way to solve the isotropic Eliashberg equations.  
In the simplest case, only the Eliashberg spectral function ``\alpha^2 F`` is required as input and it returns the superconducting critical temperature (Tc). There is also the option to save the components of the self-energy ``\Delta, Z, \chi`` at each temperature.
For more advanced calculations files containing the DOS, the screened Coulomb interaction W and the energy grid points of W are needed.
The package is able to solve the isotropic eliashberg equations within one of the following approximations:
- constant Dos + Morel-Anderson pseudopotential ``\mu^*``
- variable Dos + Morel-Anderson pseudopotential ``\mu^*``
- constant Dos + screened Coulomb interaction ``W``
- variable Dos + screened Coulomb interaction ``W``

Furthermore, results based on the Allen-Dynes and Allen-Dynes-McMillan formulas are provided.  (cite our paper??)  

## Installation
In order to run the code you need to [install](https://julialang.org/downloads/) julia 1.10 or higher.
IsoME.jl is not yet a registered package but it can be installed using the link to the github repository.
```julia-repl 
julia> using Pkg
julia> Pkg.add(path="https://github.com/cheil/IsoME.jl")
```

 ----
 **ToDo**  
IsoME.jl is a registered Julia package and can thus be installed via
```julia-repl 
julia> using Pkg
julia> Pkg.add("IsoME")
```

 ---

 ## Usage
After adding the package to your environment it can be loaded via
```julia-repl 
julia> using IsoME
```
To search for the Tc within the constant dos approximation using the Anderson-Pseuodopotential, only the path to the ``\alpha^2F``-file has to be specified. We recommend to set also the path to the output directory, beacuse per default the results will be written into the current working directory.
```julia-repl 
julia> inp = arguments(
                a2f_file="Path to a2f-file", 
                outdir="Path to output directory"
                )
```
All other input values are optional and are contained within the custom data type `arguments()`. The input arguments can be viewed by using a dot after inp, e.g.
```julia-repl 
julia> inp.omega_c
10000.0
```
gives the cutoff of the Matsubara summation in mEv.  
Finally, the calculation can be started via
```julia-repl 
julia> EliashbergSolver(inp)
```
Per default a log-file and a file containing the summary of the run, as well as figures of the superconducting gap and the a2f-values are created. 
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
Some of the input files have to obey a certain structure. For more information please refer to [Input Parameters](@ref)


## Minimal example
A minimal example can be found at [Example files](https://github.com/cheil/IsoME.jl/tree/main/test/Nb).
If you have followed the steps above you should be able to run the file test.jl.
```console
~ $ julia test.jl
```


 ## Outline 
 ```@index
```

```@contents
```

