# Best practices
## Set up a project environment
We recommend to use project [environments](https://pkgdocs.julialang.org/v1/environments/) in julia. This determines all dependencies of a project and ensures reproducibility of the results.
Furthermore, this will prevent incompatibilities with any other previously installed packages.
An environment can be set up directly in a julia REPL
```julia-repl
(@v1.9) pkg> activate MyProject
Activating new environment at `~/MyProject/Project.toml`

(MyProject) pkg> st
    Status `~/MyProject/Project.toml` (empty project)

(MyProject) pkg> add IsoME
```
Now you should have an environment containing only the IsoME package. 
When running a script test.jl from the terminal, either specify the path to the desired environment via
```console
~ $ julia --project=/MyProject/ test.jl
```
or copy the Manifest.toml and Project.toml of the environment into the folder of test.jl.
The environment can also be activated already at the beginning of the test.jl file:
```julia
    using Pkg
    
    Pkg.activate("/path/to/environment/")
```


## Do not reuse the input structure
Some parameters of the input structure may be overwritten during a run.  
Thus, we highly recommend to always set up a new instance of the input structure.  
This can be achieved by explicitly calling `arguments()` before each run. 

Lets assume you want to calculate the Tc for two different values for $\mu^*$.  
The naive approach would be to initialze the input structure once and just overwrite the $\mu^*_{AD}$ value for the second run.

!!! error "Wrong"
    ```julia
    inp = arguments(some input, muc_AD = 0.12)
    EliashbergSolver(inp)
    inp.muc_AD = 0.14
    EliashbergSolver(inp)
    ```

The output of these two runs is astonishingly exactly the same.
Reason for that is, that the value for $\mu^*_{ME}$ has been overwritten in the first run and in the second run, the code assumes that a $\mu^*_{ME}$ value has been specified manually and does not overwrite.   
To make this work, the $\mu^*_{ME}$ value has to be reseted as well.

!!! warning "Not recommended"
    ```julia
    inp = arguments(some input, muc_AD = 0.12)
    EliashbergSolver(inp)
    inp.muc_AD = 0.14
    inp.muc_ME = -1
    EliashbergSolver(inp)
    ```

The most convenient and recommended way to do this is just by overwriting the whole instance

!!! tip "Recommended"
    ```julia
    inp = arguments(some input, muc_AD = 0.12)
    EliashbergSolver(inp)
    inp = arguments(some input, muc_AD = 0.14)
    EliashbergSolver(inp)
    ```
By using this strategy, it is impossible to hand-over any unexpected input to the `EliashbergSolver()`.


## Convergence 
IsoME was designed as a robust and user-friendly framework for calculating superconducting properties. Nevertheless, for genuine Tc predictions and proper interpretation of the results, careful conducted convergence tests have to be performed.  

- **Input files:**
Accurate results can only be achieved through carefully conducted convergence tests for the input files. In particular, the ``\alpha^2 F`` data needs to be of sufficient quality. Different Brillouin zone grids or smearings can have a huge impact on the ``\mathrm{T}_C``. We highly recommend to always check the convergence for different smearings.  
- **Convergence parameters in IsoME:**
Considerable effort has been invested in selecting default parameters that, in most cases, ensure both computational efficiency and robust convergence.
Nevertheless, convergence should always be checked.
Convergence parameter include the Matsubara cutoff *omega_c* and the energy cutoff *encut*.
In both cases, the ideal cutoff is bounded from above as the adaption formula of *muc_ME* breaks down for very large *omega_c* and arbitrary high *encut*'s are against the spirit of the isotropic approximation.

Furthermore, the energy gird around the Fermi surface must be sufficiently dense. The steps and interpolation boundaries can be adapted through *itpStepSize* and *itpBounds*.


## Ab-initio calculations with ``\mu^*``
The choice of μ significantly influences the results. Traditionally, ``\mu^*`` is treated as an adjustable parameter and typically chosen within the range of 0.1 to 0.16 to fit experimental
values.  
For fully ab-initio calculations, ``\mu`` must be computed via
```math
\mu = N(\varepsilon_F)W(\varepsilon_F,\varepsilon_F)~,
```
 and ``\mu^*_E`` adapted according to 
```math
\mu^*_{E}=\frac{\mu}{1+\mu \text{ ln}\left(\frac{\varepsilon_{el}}{\hbar \omega_{c}}\right)}~.
```