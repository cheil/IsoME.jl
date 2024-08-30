# isoME
The package isoME.jl offers a quick way to solve the isotropic Eliashberg equations.  
The superconducting critical temperature (Tc) as well as the components of the self-energy $\Delta, Z, \chi$ are returned. (cite our paper??)  
In the simplest case only a file containing the Eliashberg spectral function $\alpha^2 F$ is required as input.
For more advanced calculations files containing the DOS, the screened Coulomb interaction W and the energy grid points of W are needed.
The package is able to solve the isotropic eliashberg equations within one of the following approximations:
- constant Dos + Anderson pseudopotential mu*
- full bandwidth + Anderson pseudopotential mu*
- constant Dos + screened Coulomb interaction W
- full bandwidth + screened Coulomb interaction W

Furthermore, results based on the Allen-Dynes and Allen-Dynes-McMillan formulas are provided. 

## Installation
Since isoME.jl is not yet registered, it has to be installed using the link to the github repository.

    using Pkg
    Pkg.add(path="https://github.com/cheil/IsoME.jl")

 ----
 **ToDo**  
isoME.jl is a registered Julia package and can thus be installed via

    using Pkg
    Pkg.add("isoME")
 ---

 ## Example Usage
To get started just load the package via

    using IsoME

To search for the Tc within the constant dos approximation using the Anderson-Pseuodopotential one has to specify only the path to the a2f-file. But we recommend to set also the path to the output directory, otherwise the output will be written to the current working directory.

    inp = arguments(
        a2f_file="Path to a2f-file", 
        outdir="Path to output directory"
        )

All other input values are optional and are contained within the inp-structure, where inp is of the custom type arguments. The input arguments can be viewed by using a dot after inp, e.g.

    inp.omega_c

gives the cutoff of the Matsubara summation in mEv.  
Finally, the calculation can be started via

    EliashbergSolver(inp)

Per default the results are written into an output-file, a figure of the superconducting gap as well as the a2f-values are plotted and the log-file contains further information about the run.  
Calculations within one of the other approximations can be started completely analog by specifying the corresponding inputs.
Some of the input files have to obey a certain structure. For more information please refer to LINK TO INPUT PAREMETER PAGE

### Variable DoS Anderson-Pseudopotential
    inp = arguments(
                    a2f_file            = "Path to a2f-file", 
                    outdir              = "Path to output directory",
                    dos_file            = "Path to dos-file",
                    cDOS_flag           = 0,
                    include_Weep        = 0
                )

### Constand DoS full screened Coulomb interaction
    inp = arguments(
                    a2f_file            = "Path to a2f-file", 
                    outdir              = "Path to output directory",
                    dos_file            = "Path to dos-file",
                    Weep_file           = "Path to W-file",
                    Wen_file            = "Path to W energies-file",
                    cDOS_flag           = 1,
                    include_Weep        = 1
                )

### Constand DoS full screened Coulomb interaction
    inp = arguments(
                    a2f_file            = "Path to a2f-file", 
                    outdir              = "Path to output directory",
                    dos_file            = "Path to dos-file",
                    Weep_file           = "Path to W-file",
                    Wen_file            = "Path to W energies-file",
                    cDOS_flag           = 0,
                    include_Weep        = 1
                )



 ## Outline ??

