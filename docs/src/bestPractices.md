# Best practices
## Set up a project environment
We recommend using project [environments](https://pkgdocs.julialang.org/v1/environments/) in Julia. This determines all dependencies of a project and ensures reproducibility of the results.
Furthermore, this will prevent incompatibilities with any other previously installed packages.
An environment can be set up directly in a Julia REPL
```julia-repl
(@v1.10) pkg> activate MyProject
Activating new environment at `~/MyProject/Project.toml`

(MyProject) pkg> st
    Status `~/MyProject/Project.toml` (empty project)

(MyProject) pkg> add IsoME
```
Now you should have an environment containing only the IsoME package. 
When running a script test.jl from the terminal, either specify the path to the desired environment via
```console
~ $ julia --project=~/MyProject/ test.jl
```
or copy the Manifest.toml and Project.toml of the environment into the folder of test.jl.
The environment can also be activated already at the beginning of the test.jl file:
```julia
using Pkg

Pkg.activate("/path/to/environment/")
```


## Do not reuse the input structure
Some parameters of the input structure may be overwritten during a run.  
Thus, we highly recommend always setting up a new instance of the input structure.  
This can be achieved by explicitly calling `arguments()` before each run. 

Let's assume you want to calculate the ``\mathrm{T}_C`` for two different values of ``\mu^*``.  
The naive approach would be to initialize the input structure once and just overwrite the ``\mu^*_{AD}`` value for the second run.

!!! error "Wrong"
    ```julia
    inp = arguments(some input, muc_AD = 0.12)
    EliashbergSolver(inp)
    inp.muc_AD = 0.14
    EliashbergSolver(inp)
    ```

Both of these runs will produce exactly the same Eliashberg results.
The reason is that the value of ``\mu^*_{ME}`` is overwritten during the first run. In the second run, the same ``\mu^*_{ME}`` is used as it is assumed to be user-defined.
Only the Allen-Dynes estimates differ, since those are evaluated from ``\mu^*_{AD}``.
To ensure correct behavior, the ``\mu^*_{ME}`` value must therefore be reset as well.

!!! warning "Not recommended"
    ```julia
    inp = arguments(some input, muc_AD = 0.12)
    EliashbergSolver(inp)
    inp.muc_AD = 0.14
    inp.muc_ME = NaN
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
By using this strategy, it is impossible to hand over any unexpected input to the `EliashbergSolver()`.
The same applies to `RealAxisSolver()`, which is driven by the same input structure.


## Convergence
IsoME was designed as a robust and user-friendly framework for calculating superconducting properties. Nevertheless, for genuine ``\mathrm{T}_C`` predictions and proper interpretation of the results, carefully conducted convergence tests have to be performed.  

- **Input files:**
Accurate results can only be achieved through carefully conducted convergence tests for the input files. In particular, the ``\alpha^2 F`` data needs to be of sufficient quality. Different Brillouin zone grids or smearings can have a huge impact on the ``\mathrm{T}_C``. We highly recommend always checking the convergence for different smearings.
If ``\alpha^2F`` and the DOS were computed with different ``N(\varepsilon_F)``, `a2f_Nef` puts them back on the same normalization — an expert flag that most calculations do not need.  
- **Convergence parameters in IsoME:**
Considerable effort has been invested in selecting default parameters that, in most cases, ensure both computational efficiency and robust convergence.
Nevertheless, convergence should always be checked.
On the imaginary axis the convergence parameters are the frequency cutoff `omega_c` (the Matsubara cutoff here) and the energy cutoff `encut`.
In both cases, the ideal cutoff is bounded from above as the adaptation formula of `muc_ME` breaks down for very large `omega_c` and arbitrarily high `encut` values are incompatible with the isotropic approximation.
Note that in vDOS calculations `encut` is clamped to `omega_c`, so a wider ``\varepsilon``-window requires a larger `omega_c` as well.

Furthermore, the energy grid around the Fermi level must be sufficiently dense. The steps and interpolation boundaries can be adapted through `itpStepSize` and `itpBounds`.

The real-axis solver shares `omega_c` — there it is the cutoff of the ``\omega``-grid — and adds the frequency step `domega` (which also sets the resolution of the kernel), the number of Chebyshev points per pole `n_cheb` and, in vDOS calculations, the energy step `depsilon`.
It reacts far more sensitively to these than the imaginary-axis solver does to its own — see [Real-axis grids and the kernel](@ref) and [The μ-update](@ref) on the [Troubleshooting](@ref) page for the symptoms of each.

- **Check the log of a finished run:**
A run that ends without an error is not automatically a converged one.
Search `log.txt` for `Convergence not achieved within` — a temperature that ran out of `N_it` iterations is reported there and its gap should not be trusted.
It is also worth scanning the error column of the iteration tables: it should fall steadily, and an order parameter that oscillates between positive and negative values signals too aggressive mixing (see the [FAQ](@ref)).


## Ab-initio calculations with ``\mu^*``
The choice of ``\mu`` significantly influences the results. Traditionally, ``\mu^*`` is treated as an adjustable parameter and typically chosen within the range of 0.1 to 0.16 to fit experimental values.  
For fully ab-initio calculations, ``\mu`` must be computed from
```math
\mu = N(\varepsilon_F)W(\varepsilon_F,\varepsilon_F)~,
```
 and ``\mu^*_{ME}`` adapted according to 
```math
\mu^*_{ME}=\frac{\mu}{1+\mu \ln\left(\frac{\varepsilon_{el}}{\omega_{c}}\right)}~.
```

In IsoME this happens automatically as soon as a `Weep_file` is given: ``\mu`` is evaluated from ``W`` at the Fermi level and converted to ``\mu^*_{AD}`` and ``\mu^*_{ME}``.
Setting `mu` directly has the same effect.
Either way the conversion needs a typical electronic energy ``\varepsilon_{el}`` — `typEl`, or the Fermi energy (`ef` / `efW`) when `typEl` is unset — and the run stops with an error if that combination drives ``\mu^*`` negative.
