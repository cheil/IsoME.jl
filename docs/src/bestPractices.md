# Best practices
## Set up a project environment
We recommend to use project [environments](https://docs.julialang.org/en/v1/manual/code-loading/#Environments-1) in julia. This determines all dependencies of a project and ensures reproducibility of the results. An environment can be set up via a julia REPL
```julia-repl
(@v1.9) pkg> activate MyProject
Activating new environment at `~/MyProject/Project.toml`

(MyProject) pkg> st
    Status `~/MyProject/Project.toml` (empty project)

(MyProject) pkg> add IsoME
```
Now you should have an environment containing only the IsoME package. 
When running a script test.jl from the terminal, it can be necessary to specify the environment you want to use explicitly
```console
~ $ julia --project=/MyProject/ test.jl
```


## Do not reuse the input structure
Some parameters of the input structure may be overwritten during a run.  
Thus, we highly recommend to always set up a new instance of the input structure.  
This can be achieved by just explicitly calling `arguments()` before each run. 

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
Reason for that is, that the value for $\mu^*_{ME}$ has been overwritten in the first run. So in the second run, the code assumes that you have entered a $\mu^*_{ME}$ value by hand and does not change it.   
To make this work, you would have to reset the $\mu^*_{ME}$ value as well.

!!! warning "Not recommended"
    ```julia
    inp = arguments(some input, muc_AD = 0.12)
    EliashbergSolver(inp)
    inp.muc_AD = 0.14
    inp.muc_ME = -1
    EliashbergSolver(inp)
    ```

But the most convenient way to do this is just by overwriting the whole instance

!!! tip "Recommended"
    ```julia
    inp = arguments(some input, muc_AD = 0.12)
    EliashbergSolver(inp)
    inp = arguments(some input, muc_AD = 0.14)
    EliashbergSolver(inp)
    ```
By using this strategy, it is impossible to hand-over an unexpected input to the `EliashbergSolver()`.


