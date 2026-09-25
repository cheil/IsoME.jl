# Input
Input parameters are collected in the [composite type](https://docs.julialang.org/en/v1/manual/types/#Composite-Types) `arguments`.
Only the path to the ``\alpha^2F`` file is mandatory. Everything else either has a default or is inferred during the run.
Fields that are meant to be inferred are left at a sentinel — `NaN` (real-valued fields), `-1` (integer fields such as line counts and column indices) or `""` (strings) — and are overwritten once their value is known.
Before IsoME 2.0 the real-valued fields used `-1` for this, so passing `-1` to one of them is now rejected with a message pointing at `NaN`.
`mu`, `muc_AD`, `muc_ME`, `mixing_beta`, `conv_thr`, `minGap`, `N_it` and `min_it` are likewise rejected if given a negative value, both when the struct is built and again at the start of a solve. `Tc_tol` must be strictly positive, and `N_it` must exceed both `min_it` and `nItFullCoul + 1`, since convergence is only accepted after both have passed.
Because of this, we recommend creating a fresh `arguments` instance for each call to `EliashbergSolver()` or `RealAxisSolver()` — see [Best practices](@ref).

All energies are handled internally in meV.
If an input file uses another supported unit, IsoME tries to extract the unit from the file header and convert the data automatically.
Should the automatic detection fail, check the header of the input file or set the corresponding unit parameter manually.

The inputs are grouped below into general inputs, inputs specific to the imaginary-axis solver, and inputs specific to the real-axis solver.
The last two groups are described in detail on the [Matsubara Solver](@ref) and [Real Axis Solver](@ref) pages.
The formats of the input files themselves are documented in the [Read-In](@ref) section further down.

## General inputs
These inputs are shared by both solvers unless noted otherwise.

| Name    |      Type      |   Default   | Description | Comment  |
|:--------|:---------------|:------------|:------------|:---------|
| temps   | Vector{Float64} |  [-1.0] | Temperatures considered in the calculation | `[-1.0]`: search for ``\mathrm{T}_C``; otherwise solve at the specified temperatures |
| a2f_file    | String   |  -  | Path to the ``\alpha^2F`` file | The only mandatory input |
| ind_smear   | Int64 |  -1 | Smearing column used from the ``\alpha^2F`` file | `-1`: the middle column is used |
| cDOS_flag | Int64 |  1   | Selects constant or variable DOS mode | 0: variable DOS; 1: constant DOS |
| include_Weep | Int64 | 0 | Selects Morel-Anderson or screened-Coulomb mode | 0: ``\mu^*`` approximation; 1: static ``W(\varepsilon,\varepsilon')`` interaction |
| dos_file  |  String  |     ""    | Path to the DOS file | Required if `cDOS_flag = 0` or `include_Weep = 1` |
| Weep_file |  String  |     ""    | Path to the ``W`` file | Required if `include_Weep = 1` |
| Wen_file  |  String  |     ""    | Path to a file containing the energy grid points of ``W`` | Only required if the grid is not contained in `Weep_file` |
| mu      | Float64 |   NaN  | ``\mu=N(\varepsilon_F)W(\varepsilon_F,\varepsilon_F)`` | Measure of the Coulomb strength; `NaN`: inferred |
| muc_AD  | Float64 |  NaN | Morel-Anderson pseudopotential for Allen-Dynes estimates, ``\mu^*_{AD}`` | `NaN`: the default convention described below is used |
| muc_ME  | Float64 | NaN | Morel-Anderson pseudopotential for Migdal-Eliashberg calculations, ``\mu^*_{ME}`` | `NaN`: inferred from `muc_AD`, `mu` or ``W`` when possible |
| typEl | Float64 |  NaN | Typical electronic energy scale | In meV; used to calculate ``\mu^*`` from ``\mu`` |
| ef        | Float64 |  NaN | Fermi energy of the DOS | In meV; `NaN`: extracted from the DOS-file header |
| efW       | Float64 |  NaN | Fermi energy of the ``W`` grid | In meV; `NaN`: extracted from the ``W``-file header |
| omega_c | Float64 | 7000.0 | Frequency cutoff ``\omega_c`` | In meV; the Matsubara cutoff on the imaginary axis and the ``\omega``-grid cutoff on the real axis |
| encut  | Float64 |  2000.0  | Symmetric energy cutoff of the ``\varepsilon``-grid | In meV; the grid spans ``\pm`` `encut` and bounds every ``\varepsilon``-integration, including ``\chi`` and the charge neutrality condition. In vDOS runs it is clamped to `omega_c`, on the imaginary and the real axis alike |
| mu_flag    | Int64 | 1   | Update the chemical potential in vDOS calculations | 0: no; 1: yes, recommended |
| mixing_beta | Float64 | NaN | Linear mixing factor | `NaN`: use the default iteration-dependent schedule. Both solvers mix linearly |
| nItFullCoul | Int64 |  10    | Number of iterations used to ramp up the Coulomb contribution | Helps stabilize the initial iterations |
| conv_thr | Float64 | ``10^{-4}`` | Convergence threshold | Applied to the gap update |
| minGap   | Float64 | 0.1  | Lower gap threshold | In meV; the temperature is treated as normal conducting if the gap drops below `minGap` |
| N_it | Int64 | 5000 | Maximum number of iterations | - |
| min_it | Int64 | 10 | Minimum number of iterations before convergence is accepted | Used by both solvers; the Coulomb ramp (`nItFullCoul`) must have finished as well |
| Tc_tol | Float64 | 1.0 | Resolution the ``\mathrm{T}_C`` search is carried to | In K; must be strictly positive. Once the search brackets ``\mathrm{T}_C``, it refines on a lattice of this spacing until the bracket is `Tc_tol` wide, so the reported ``\mathrm{T}_C`` carries an uncertainty of ``\pm`` `Tc_tol`/2. Smaller values cost additional temperatures; ignored when `temps` is given explicitly. `Tc_tol` only narrows the bracket: a temperature counts as normal conducting once the gap drops below `minGap` or the iteration fails to converge within `N_it`, so the bracket closes around the temperature where that happens, which lies somewhat below the true ``\mathrm{T}_C``. Refining far below that offset does not make ``\mathrm{T}_C`` more accurate |
| outdir | String | `joinpath(pwd(), "IsoME")` | Path to the output directory | An existing directory is never written into: a run counter is appended, so repeated runs land in `IsoME_1/`, `IsoME_2/`, and so on |
| flag_figure | Int64 |  1 | Plot the gap and ``\alpha^2F`` values | 0: no; 1: yes |
| flag_writeSelfEnergy | Int64 | 0  | Save **and** plot the self-energy components (`.dat` + `.png`) | 0: no; 1: yes. Works for both solvers and all modes; files are written to `outdir/SelfEnergy/`. The first header line of each `.dat` carries the converged chemical potential as `mu_F = … meV` (relative to the ``\varepsilon_F`` of the DOS input), which any ``\varepsilon``-resolved post-processing needs |
| material | String | "Material" | Name of the compound | Used in plots and summaries |
| returnTc | Bool | false | Return the estimated ``\mathrm{T}_C`` interval | Mainly useful for scripts and tests |
| testMode | Bool | false | Suppress all file output | Used by the test suite |


## Imaginary-axis inputs
These inputs control the imaginary-axis solver `EliashbergSolver()` and its post-processing.

| Name    |      Type      |   Default   | Description | Comment  |
|:--------|:---------------|:------------|:------------|:---------|
| sparseSamplingTemp | Float64 | 2.0 | Temperature below which sparse sampling is used | In K |
| itpBounds | Vector{Float64} | [100.0, 500.0] | Bounds of the DOS interpolation regions around the Fermi level | In meV |
| itpStepSize | Vector{Int64} | [1, 5, 50] | Step sizes used within the interpolation regions | In meV; one entry more than `itpBounds` |
| flag_acon | Bool | false | Analytically continue the converged solution to real frequencies and plot it | Padé approximants; imaginary axis only. Written to `outdir/ACON/`; the real-frequency window is ``\pm`` `omega_c` with a step of 1 meV. Independent of `flag_writeSelfEnergy` |

## Real-axis inputs
These inputs control the direct real-axis solver `RealAxisSolver()`, which supports cDOS``+\mu``, vDOS``+\mu`` and vDOS``+W``.
Note that `include_Weep = 1` requires `cDOS_flag = 0` here.
The ``\varepsilon``-grid of the real-axis solver is uniform (step `depsilon`) and does not use `itpBounds`/`itpStepSize`.
The cutoff of the ``Z(\omega)``, ``\Delta(\omega)`` and ``\chi(\omega)`` grids — which also bounds the ``\omega'``-integration — is the shared `omega_c` listed under the general inputs; in vDOS, `encut` is clamped to at most that value, on this axis and on the imaginary one alike.

| Name    |      Type      |   Default   | Description | Comment  |
|:--------|:---------------|:------------|:------------|:---------|
| domega | Float64 | 1.0 | Step of the ``Z(\omega)``, ``\Delta(\omega)`` and ``\chi(\omega)`` grids | In meV; also the step of the ``\omega'``-tail **and** of the tabulated kernel |
| n_cheb | Int64 | 1000 | Chebyshev points per pole in the head region of the ``\omega'``-grid | Resolves the poles of the ``\omega'``-integrand |
| depsilon | Float64 | 10.0 | Step of the ``\varepsilon``-grid | In meV; vDOS calculations only |

The ``\omega``-grid does not start at ``0``: its first point is the smallest multiple of `domega` that is at least ``0.1`` meV, because several integrands behave like ``1/\omega``.
`domega` is a single step size that all real-axis grids share — the ``\omega``-grids, the ``\omega'``-tail and the kernel tables — so it is the one convergence parameter of the frequency discretization.
The kernel tables are only ever sampled at multiples of `domega`, which is why they carry no step of their own.


### Pseudopotentials ``\mu,~\mu^*_{AD}~\&~\mu^*_{ME}``
``\mu`` measures the strength of the Coulomb interaction at the Fermi surface: ``\mu=N(\varepsilon_F)W(\varepsilon_F,\varepsilon_F)``.
It is connected to the pseudopotentials via:
```math
\mu^*_{AD}=\frac{\mu}{1+\mu \ln\left(\frac{\varepsilon_{el}}{\omega_{ph}}\right)}
```
where ``\omega_{ph}`` is a characteristic cutoff frequency for the phonon-induced interaction and ``\varepsilon_{el}`` is a characteristic electronic energy scale.
The typical electronic energy can be specified explicitly through `typEl`; otherwise, the Fermi energy (`ef` or `efW`) will be used.
For the characteristic phonon cutoff, AD uses the largest frequency at which ``\alpha^2F`` rises above ``10^{-2}`` (not the largest frequency in the file), while ME uses the frequency cutoff `omega_c`, which is shared by both solvers.

By default, the ``\mu`` and ``\mu^*`` values are unset (`NaN`), which means that `muc_AD` = 0.12 is used. This also fixes `muc_ME` through
```math
\mu^*_{ME}= \frac{\mu^*_{AD}}{(1 + \mu^*_{AD} \ln(\frac{\omega_{ph}}{\omega_c}))}~.
```
If ``\mu`` or one of the ``\mu^*`` values is set, the remaining ``\mu^*`` values are calculated from it. However, no user input is overwritten.
Furthermore, if neither ``\mu`` nor ``\mu^*`` is specified but a `Weep_file` is given, ``\mu`` is calculated from ``W``.


### Interpolation of the energy grid
Since the electronic structure close to the Fermi level drives the results, the imaginary-axis solver interpolates the DOS (and ``W``) onto a piecewise-uniform energy grid defined by `itpBounds` and `itpStepSize`.
The bounds define regions around the Fermi level, and each region uses the corresponding step size from `itpStepSize`; the first step size is used within the first region around the Fermi level, the last one beyond the outermost bound.
With the defaults, the grid uses a 1 meV step within ``\pm 100`` meV of ``\varepsilon_F``, 5 meV up to ``\pm 500`` meV, and 50 meV out to `encut`.

The real-axis solver instead uses a uniform grid of step `depsilon` over the whole `encut` range, because its ``\varepsilon``-integrals are evaluated analytically on a piecewise-linear DOS.


## Read-In
IsoME automatically recognizes the format of QE, EPW and BerkeleyGW files.
Compatibility with other DFT/DFPT/GW packages is currently under development.
However, the read-in function has been designed to rely on as little formatting as possible and will often work for other formats as well.
If the auto-recognition fails, either adapt the format of the input files or set it manually through the dedicated flags.
For details please refer to the dedicated section for each input file.
Possible sources of errors are the number of header/footer lines, the Fermi energy, or the units.

In its auto-recognition mode, IsoME interprets non-numeric rows at the beginning and end of the file as the header and footer, respectively. From the header, the unit (currently supported: meV, eV, THz, Ry, Ha) and, for DOS or ``W`` files, the Fermi energy are extracted.


### ``\alpha^2F``
The Eliashberg spectral function ``\alpha^2F(\omega)`` is required for all calculations.
The ``\alpha^2F(\omega)`` -file can contain an arbitrary number of columns, but the first column must contain the energies and the remaining columns are interpreted as ``\alpha^2F(\omega)`` -values for different smearings. If the user does not specify the smearing column via `ind_smear`, the column in the middle will be used.


#### Summary formatting:
- **header:** Non-numeric rows at the beginning of the document. If the header contains the unit (meV, eV, THz, Ry, Ha), it is extracted automatically; otherwise, set the unit via `a2f_unit`.
- **footer:** Non-numeric rows at the end of the document.
- **first column:** energies
- **second column onwards:** ``\alpha^2F`` values for different smearings. By default, the smearing in the middle is used.

The number of header/footer lines and smearing values should be recognized automatically. If this is not the case, set them through the dedicated input parameters:

|     Name     |  Type  |  Default  |     Description       |     Comment         |
|--------------|--------|:---------:|-------------------------------------|---------------------------|
| a2f_unit     | String |    ""     | Unit of the a2f-file                | Extracted from the header of the a2f-file if not set. Currently supported: meV, eV, THz, Ry, Ha |
| nsmear       | Int64  |    -1     | number of smearings in the a2f-file | Auto-recognition if unset |
| nheader_a2f  | Int64  |    -1     | number of header lines in a2f_file  | Auto-recognition if unset |
| nfooter_a2f  | Int64  |    -1     | number of footer lines in a2f_file  | Auto-recognition if unset |
| a2f_Nef      | Float64|    NaN    | ``N(\varepsilon_F)`` used when ``\alpha^2F`` was computed | Optional. If set, ``\alpha^2F`` is rescaled to the DOS-file ``N(\varepsilon_F)`` (see below) |

!!! details "Detailed description"
    >  **a2f_file** :: STRING
    >
    > Path to the ``\alpha^2F``-file.
    > Required for all calculations.
    > The first column must contain the energies, the second column onwards the ``\alpha^2F`` values for different smearings.
    >

    > **ind_smear** :: INTEGER | Default: -1
    >
    > Index of the smearing that should be used. If unset, the column in the middle is used.
    >

    >  **a2f_unit** :: STRING | Default: ""
    >
    > Energy unit in the ``\alpha^2F``-file
    >
    > | Value | Description |
    > | ---   | ----------- |
    > | ""   | Auto-extraction from header |
    > | meV     | - |
    > | eV  | - |
    > | THz     | - |
    > | Ry     | - |
    > | Ha     | - |

    >  **nheader_a2f** :: INTEGER | Default: -1
    >
    > Number of header lines in the ``\alpha^2F``-file
    >
    > | Value | Description |
    > | ---   | ----------- |
    > | -1 | Auto-recognition |

    >  **nfooter_a2f** :: INTEGER  | Default: -1
    >
    > Number of footer lines in the ``\alpha^2F``-file
    >
    > | Value | Description |
    > | ---   | ----------- |
    > | -1 | Auto-recognition |

    >  **a2f_Nef** :: FLOAT | Default: NaN
    >
    > Density of states at the Fermi level, ``N(\varepsilon_F)``, that was used when the
    > ``\alpha^2F``-file was computed. Optional; leave at `NaN` to disable.
    >
    > When set, ``\alpha^2F`` is rescaled by ``a2f\_Nef / N(\varepsilon_F)``, where
    > ``N(\varepsilon_F)`` is read from the DOS file. This puts ``\alpha^2F`` and the electronic
    > normalization used in the Eliashberg equations on the same ``N(\varepsilon_F)``, and it is
    > applied once during read-in so that every downstream quantity (``\mu^*``, the Allen-Dynes
    > estimate, and the solver) uses the rescaled ``\alpha^2F``. The applied factor is printed to the
    > console and log.
    >
    > A DOS file is required: in vDOS (and Weep+``\mu``) runs the value already read is used; in
    > cDOS+``\mu`` runs ``N(\varepsilon_F)`` is read from `dos_file` if one is given. If no DOS file is
    > available the rescaling is skipped with a warning.
    >
    > **Unit convention:** `a2f_Nef` must be given in the same convention the DOS file is reduced to
    > internally: per meV and per spin, i.e. after the division by `spinDos`. Check the printed factor:
    > it should be of order 1.


### ``N(\epsilon)``
For vDOS calculations a DOS file is required.
The first and second columns of the DOS file are interpreted as the energies and DOS values, respectively. All other columns are ignored. The DOS values are divided by `spinDos` (2 by default) to remove the double counting from spin degeneracy. If this is not desired, set `spinDos` to 1.
The energy grid does not have to be uniform: the DOS is interpolated linearly on the energies as given. The energies must be in ascending order without repetitions; otherwise the read-in stops with an error.

#### Summary formatting:
- **header:** Non-numeric rows at the beginning of the document. If the header contains the unit (meV, eV, THz, Ry, Ha), it is extracted automatically; otherwise, set the unit via `dos_unit`. If the header contains only one numeric value, this is interpreted as the Fermi energy. If there are several numerical values, IsoME checks for a keyword (`ef`, `efermi`, ...) indicating the Fermi energy. If extraction fails, adapt the header or set the Fermi energy via `ef`.
- **footer:** Non-numeric rows at the end of the document.
- **first column:** energies
- **second column:** dos values

|     Name     |  Type  |  Default  |          Description          |                Comment              |
|--------------|--------|:---------:|---------------------------|-----------------------------------------|
| spinDos      | Int64  |     2     | Whether the DOS includes spin degeneracy | 1 = spin not considered <br> 2 = spin considered |
| dos_unit     | String |    ""     | Unit in the DOS file | Extracted from the header if not set <br> Currently supported: meV, eV, THz, Ry, Ha |
| nheader_dos  | Int64  |    -1     | Number of header lines in the DOS file | Auto-recognition if unset |
| nfooter_dos  | Int64  |    -1     | Number of footer lines in the DOS file | Auto-recognition if unset |

!!! details "Detailed description"
    >  **dos_file** :: STRING
    >
    > Path to the DOS file
    > Only required for vDOS calculations

    >  **spinDos** :: INTEGER | Default: 2
    >
    > Spin convention used in the DOS file
    >
    > | Value | Description |
    > | ---   | ----------- |
    > | 1     | Spin is not considered in the DOS |
    > | 2     | DOS contains spin degeneracy |

    >  **dos_unit** :: STRING | Default: ""
    >
    > Energy unit in the DOS file
    >
    > | Value | Description |
    > | ---   | ----------- |
    > | ""   | Auto-extraction from header |
    > | meV     | - |
    > | eV     | - |
    > | THz     | - |
    > | Ry     | - |
    > | Ha     | - |

    >  **nheader_dos** :: INTEGER | Default: -1
    >
    > Number of header lines in the DOS file
    >
    > | Value | Description |
    > | ---   | ----------- |
    > | -1 | Auto-recognition |

    >  **nfooter_dos** :: INTEGER  | Default: -1
    >
    > Number of footer lines in the DOS file
    >
    > | Value | Description |
    > | ---   | ----------- |
    > | -1 | Auto-recognition |


### Weep
The `Weep_file` is required for ``W`` calculations.
By default, IsoME assumes that the third column contains the ``W(\varepsilon,\varepsilon')`` values and that the first and second columns contain the corresponding energy-grid coordinates.
The columns can be changed via `Weep_col` and `Wen_col`.
If the row and column identifiers are consecutive indices rather than energies, provide an additional `Wen_file` containing the ``W`` energy-grid points.
As for the DOS, the ``W`` energy grid does not have to be uniform; ``W`` is interpolated bilinearly on the grid as given, which must be in ascending order without repetitions.
#### Summary formatting:
- **header:** Non-numeric rows at the beginning of the document. If the header contains the unit (meV, eV, THz, Ry, Ha), it is extracted automatically; otherwise, set the unit via `Weep_unit`. If the header contains only one numeric value, it is interpreted as the Fermi energy. If there are several numerical values, IsoME checks for a keyword (`ef`, `efermi`, ...) indicating the Fermi energy. If extraction fails, adapt the header or set the Fermi energy via `efW`.
- **footer:** Non-numeric rows at the end of the document.
- **columns:** energy-grid points are assumed to be in the first and second columns, and the first column is used by default (`Wen_col` = 1). The ``W`` values are assumed to be in the third column by default (`Weep_col` = 3).


If the `Weep_file` does not contain the energies, an additional `Wen_file` can be specified.
- **header:** Non-numeric rows at the beginning of the document. It is assumed that the header contains the unit. If not, the unit has to be specified via `Wen_unit`.
- **footer:** Non-numeric rows at the end of the document.
- **first column:** energy grid points ``\epsilon`` of ``W(\epsilon, \epsilon')``. The column containing the energies can be changed via `Wen_col`.


|     Name     |  Type  |  Default  |          Description          |                Comment              |
|--------------|--------|:---------:|---------------------------|-----------------------------------------|
| Weep_unit    | String |     ""    | Unit of ``W``             | Extracted from the header of `Weep_file` if not set <br> Currently supported: meV, eV, THz, Ry, Ha |
| Wen_unit     | String |     ""    | Unit of the ``W`` energies    | Extracted from the header of `Wen_file` if not set <br> Currently supported: meV, eV, THz, Ry, Ha |
| Weep_col     | Int64  |     3     | Column containing the ``W`` values in `Weep_file` | |
| Wen_col      | Int64  |     1     | Column containing the ``W`` energy grid in `Weep_file` or `Wen_file` | |
| nheader_Weep | Int64  |    -1     | Number of header lines in `Weep_file` | Auto-recognition if unset |
| nfooter_Weep | Int64  |    -1     | Number of footer lines in `Weep_file` | Auto-recognition if unset |
| nheader_Wen  | Int64  |    -1     | Number of header lines in `Wen_file` | Auto-recognition if unset |
| nfooter_Wen  | Int64  |    -1     | Number of footer lines in `Wen_file` | Auto-recognition if unset |

!!! details "Detailed description"
    >  **Weep_file** :: STRING
    >
    > Path to the file containing ``W(\varepsilon,\varepsilon')``.
    > Required for ``W`` calculations.

    >  **Wen_file** :: STRING
    >
    > Optional path to a file containing the energy grid of ``W(\varepsilon,\varepsilon')``.
    > Only required if the grid is not contained in `Weep_file`.

    >  **Weep_unit** :: STRING | Default: ""
    >
    > Energy unit in the ``W`` file.
    >
    > | Value | Description |
    > | ---   | ----------- |
    > | ""   | Auto-extraction from header |
    > | meV     | - |
    > | eV      | - |
    > | THz     | - |
    > | Ry      | - |
    > | Ha      | - |

    >  **Wen_unit** :: STRING | Default: ""
    >
    > Energy unit in the optional ``W`` energy-grid file.
    >
    > | Value | Description |
    > | ---   | ----------- |
    > | ""   | Auto-extraction from header |
    > | meV     | - |
    > | eV      | - |
    > | THz     | - |
    > | Ry      | - |
    > | Ha      | - |

    >  **Weep_col** :: INTEGER | Default: 3
    >
    > Column containing the ``W`` values in `Weep_file`.

    >  **Wen_col** :: INTEGER | Default: 1
    >
    > Column containing the ``W`` energy-grid values in `Weep_file` or `Wen_file`.



## Docstring
The same overview is available from the REPL through `?arguments`.

```@docs
arguments
```


## Version
Julia 1.10 or higher is required.
