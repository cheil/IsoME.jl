# Input
Input parameters are collected in the [composite type](https://docs.julialang.org/en/v1/manual/types/#Composite-Types) `arguments()`.
We recommend creating a fresh `arguments()` instance for each call to `EliashbergSolver()` or `RealAxisSolver()`, because some fields are inferred or updated during a run.
Fields that are meant to be inferred are usually marked by `-1` for numeric values or `""` for strings.

All energies are handled internally in meV.
If an input file uses another supported unit, IsoME tries to extract the unit from the file header and convert the data automatically.
If automatic unit detection fails, check the header of the input file or set the corresponding unit parameter manually.

The most relevant inputs are grouped below into general inputs, imaginary-axis inputs, and real-axis inputs.
More detailed descriptions of the input-file formats are given in the following sections.

## General inputs
These inputs are shared by the imaginary-axis and real-axis solvers unless noted otherwise.

| Name    |      Type      |   Default   | Description | Comment  |
|:--------|:---------------|:------------|:------------|:---------|
| temps   | Vector{Number} |  [-1] | Temperatures considered in the calculation | `[-1]`: search for ``T_c``; otherwise solve at the specified temperatures |
| a2f_file    | String   |    ""    | Path to ``\alpha^2F``-file | Required for all calculations |
| ind_smear   | Int64    |    -1    | Smearing column used from the ``\alpha^2F`` file | The middle column is used by default |
| cDOS_flag | Int64 |  1   | Selects constant or variable DOS mode | 0: variable DOS; 1: constant DOS |
| include_Weep | Int64 | 0 | Selects Morel-Anderson or screened-Coulomb mode | 0: ``\mu^*`` approximation; 1: static ``W(\varepsilon,\varepsilon')`` interaction |
| dos_file  |  String  |     ""    | Path to the DOS file | Required if `cDOS_flag = 0` or `include_Weep = 1` |
| Weep_file |  String  |     ""    | Path to the ``W`` file | Required if `include_Weep = 1`; real-axis calculations currently require `include_Weep = 0` |
| Wen_file  |  String  |     ""    | Path to a file containing the energy grid points of ``W`` | Only required if the grid is not contained in `Weep_file` |
| mu      | Float64        |   -1  | ``\mu=N(\varepsilon_F)W(\varepsilon_F,\varepsilon_F)`` | Measure of the Coulomb strength |
| muc_AD     | Float64     |  -1 | Morel-Anderson pseudopotential for Allen-Dynes estimates, ``\mu^*_{AD}`` | If unset, the default convention described below is used |
| muc_ME  | Float64        | -1 | Morel-Anderson pseudopotential for Migdal-Eliashberg calculations, ``\mu^*_{ME}`` | Inferred from `muc_AD`, `mu`, or ``W`` when possible |
| typEl | Float64|  -1 | Typical electronic energy scale | In meV; used to calculate ``\mu^*`` from ``\mu`` |
| ef        | Float64  |     -1    | Fermi energy of the DOS in meV | Extracted from the DOS-file header if not set |
| efW       | Float64  |     -1    | Fermi energy of the ``W`` grid in meV | Extracted from the ``W``-file header if not set |
| mu_flag    | Int64 | 1   | Update the chemical potential in vDOS calculations | 0: no; 1: yes, recommended |
| mixing_beta | Number     | -1 | Linear mixing factor | `-1`: use the default iteration-dependent schedule |
| nItFullCoul | Number     |  10    | Number of iterations used to ramp up the Coulomb contribution | Helps stabilize the initial iterations |
| conv_thr | Float64 | ``10^{-4}`` | Convergence threshold | Applied to the gap update |
| minGap   | Float64 | 0.1  | Lower gap threshold | Stop if ``\Delta(0) <`` `minGap` |
| N_it | Int64 | 5000 | Maximum number of iterations | - |
| min_it | Int64 | 10 | Minimum number of iterations before convergence is accepted | - |
| outdir | String |  pwd() | Path to the output directory | |
| flag_figure | Int64 |  1 | Plot the gap and ``\alpha^2F`` values | 0: no; 1: yes |
| flag_writeSelfEnergy | Int64 | 0  | Save self-energy components | 0: no; 1: yes |
| material | String | "Material" | Name of the compound | Used in plots and summaries |
| returnTc | Bool | false | Return the estimated ``T_c`` interval from `EliashbergSolver()` | Mainly useful for scripts and tests |


## Imaginary-axis inputs
These inputs control the imaginary-axis solver `EliashbergSolver()` and related post-processing.

| Name    |      Type      |   Default   | Description | Comment  |
|:--------|:---------------|:------------|:------------|:---------|
| imOmega_c | Float64 | 7000.0 | Matsubara cutoff ``\omega_c`` | In meV |
| encut  | Float64 |  5000.0  | Energy cutoff for DOS and ``W`` integrations | In meV |
| shiftcut  | Float64 | 2000.0 | Energy cutoff for the shift and charge-neutrality integrations | In meV; should be smaller than `encut` |
| sparseSamplingTemp | Float64 | 2.0 | Maximum temperature for sparse sampling | Used in ``T_c`` search mode |
| itpBounds | Vector{Float64} | [100, 500] | Bounds of DOS interpolation regions around the Fermi level | In meV |
| itpStepSize | Vector{Int64} | [1, 5, 50] | Step sizes used in the interpolation regions | In meV |
| flag_acon | Bool | false | Run analytic continuation after the imaginary-axis solution | Uses the implemented real-frequency continuation routines |
| plot_flag | Bool | false | Plot additional diagnostic quantities | Intended for development and diagnostics |

## Real-axis inputs
These inputs control the direct real-axis solver `RealAxisSolver()`.
The real-axis solver currently supports the ``\mu`` approximation; set `include_Weep = 0`.
For vDOS real-axis calculations, provide a DOS file and set `cDOS_flag = 0`.

| Name    |      Type      |   Default   | Description | Comment  |
|:--------|:---------------|:------------|:------------|:---------|
| reOmega_c  | Float64 | 2000.0 | Real-frequency cutoff for ``Z(\omega)`` and ``\Delta(\omega)`` | In meV |
| numReal_c | Int64 | 5000 | Number of real-frequency grid points for ``Z(\omega)`` and ``\Delta(\omega)`` | - |
| reOmega_c_shift | Float64 | 15000.0 | Real-frequency cutoff for ``\chi(\omega)`` in vDOS calculations | In meV |
| numReal_c_shift | Int64 | 25000 | Number of real-frequency grid points for ``\chi(\omega)`` | - |
| n_cheb | Int64 | 5000 | Number of Chebyshev points used for the ``\Omega`` integration | Controls the dynamical kernel integration |
| num_wp1 | Int64 | 1000 | Number of inner ``\omega'`` grid points near the gap edge | Used in vDOS real-axis calculations |
| num_wp2 | Int64 | 5000 | Number of outer ``\omega'`` grid points | Used in vDOS real-axis calculations |
| wp_max | Float64 | 2.0 | Width/location parameter for the outer ``\omega'`` Chebyshev grid | In meV |


### Pseudopotentials ``\mu,~\mu^*_{AD}~\&~\mu^*_{ME}``
``\mu`` measures the strength of the Coulomb interaction at the Fermi surface: ``\mu=W(\varepsilon_F,\varepsilon_F)``.
It is connected to the pseudopotentials via:
```math
\mu^*_{AD}=\frac{\mu}{1+\mu \text{ ln}\left(\frac{\varepsilon_{el}}{\hbar \omega_{ph}}\right)}
```
where ``\omega_{ph}`` is a characteristic cutoff frequency for the phonon-induced interaction and ``\varepsilon_{el}`` is a characteristic electronic energy scale.
The typical electronic energy can be specified explicitly through *typEl*; otherwise, the Fermi energy (*ef* or *efW*) will be used. For the characteristic phonon cutoff, the Matsubara cutoff and the maximum given phonon frequency are used for ME and AD, respectively.

By default, the ``\mu`` and ``\mu^*`` values are set to -1, which indicates that *muc_AD* = 0.12 is used. This also fixes *muc_ME* through
```math
\mu^*_{ME}= \frac{\mu^*_{AD}}{(1 + \mu^*_{AD} \ln(\frac{\omega_{ph}}{\omega_c}))}~.
```
If ``\mu`` or one of the ``\mu^*`` values is set, the remaining ``\mu^*`` values are calculated from it. However, no user input is overwritten.
Furthermore, if neither ``\mu`` nor ``\mu^*`` is specified but a `Weep_file` is given, ``\mu`` is calculated from ``W``.


### Interpolation of the energy grid
To resolve the properties of the electron-phonon coupling, IsoME interpolates the energy grid using the bounds and steps specified by *itpBounds* and *itpStepSize*. The bounds define regions around the Fermi level, and each region uses the corresponding step size from *itpStepSize*. For example, the first step size is used within the first region around the Fermi level.


## Read-In
IsoME automatically recognizes the formatting of ``\textsc{QE/EPW/BerkeleyGW}`` files.
Compatibility with other DFT/DFPT/GW packages is currently under development. However, the read-in function has been designed to rely on as little formatting as possible and will often work for different formattings as well. If the auto-recognition fails, either adapt the formatting of the input files or manually set it through the dedicated flags. For details please refer to the dedicated section for each input file.
Possible sources of errors are the number of header/footer lines, the Fermi energy, or the units.

In its auto-recognition mode, IsoME interprets non-numerical rows at the beginning and end of the file as the header and footer, respectively. From the header, the unit (currently supported: meV, eV, THz, Ry, Ha) and, for DOS or ``W`` files, the Fermi energy are extracted.


### ``\alpha^2F``
The Eliashberg spectral function ``\alpha^2F(\omega)`` is required for all calculations.
The ``\alpha^2F(\omega)`` -file can contain an arbitrary amount of columns, but the first column must contain the energies and the remaining columns are interpreted as ``\alpha^2F(\omega)`` -values for different smearings. If the user does not specify the smearing column via *ind_smear*, the column in the middle will be used.


#### Summary formatting:
- **header:** Non-numeric rows at the beginning of the document. If the header contains the unit (meV, eV, THz, Ry, Ha), it is extracted automatically; otherwise, set the unit via *a2f_unit*.
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

!!! details "Detailed description"
    >  **a2f_file** :: STRING
    >
    > Path to the ``\alpha^2F``-file.
    > Required for all calculations.
    > The first column must contain the energies, the second column onwards the ``\alpha^2F`` values for different smearings.
    >

    > **ind_smear** :: INTEGER | Default: 1
    >
    > Index of the smearing that should be used.
    >

    >  **a2f_unit** :: STRING    | Default: nothing
    >
    > Energy unit in the ``\alpha^2F``-file
    >
    > | Value | Description |
    > | ---   | ----------- |
    > | -1   | Auto-extraction from header |
    > | meV     | - |
    > | eV  | - |
    > | THz     | - |
    > | Ry     | - |
    > | Ha     | - |

    >  **nheader_a2f** :: INTEGER | Default: nothing
    >
    > Number of header lines in the ``\alpha^2F``-file
    >
    > | Value | Description |
    > | ---   | ----------- |
    > | -1   | Auto-recognition |

    >  **nfooter_a2f** :: INTEGER  | Default: nothing
    >
    > Number of footer lines in the ``\alpha^2F``-file
    >
    > | Value | Description |
    > | ---   | ----------- |
    > | -1   | Auto-recognition |


### ``N(\epsilon)``
For vDOS calculations a DOS file is required.
The first and second columns of the DOS file are interpreted as the energies and DOS values, respectively. All other columns are ignored. The DOS values are divided by two to remove double counting due to spin. If this is not desired, set the *spinDos* flag to 1.

#### Summary formatting:
- **header:** Non-numeric rows at the beginning of the document. If the header contains the unit (meV, eV, THz, Ry, Ha), it is extracted automatically; otherwise, set the unit via *dos_unit*. If the header contains only one numeric value, this is interpreted as the Fermi energy. If there are several numerical values, IsoME checks for a keyword (`ef`, `efermi`, ...) indicating the Fermi energy. If extraction fails, adapt the header or set the Fermi energy via *ef*.
- **footer:** Non-numeric rows at the end of the document.
- **first column:** energies
- **second column:** dos values

|     Name     |  Type  |  Default  |          Description          |                Comment              |
|--------------|--------|:---------:|---------------------------|-----------------------------------------|
| spinDos      | Int64  |     2     | Whether the DOS includes spin degeneracy | 1 = spin not considered ``\\`` 2 = spin considered |
| dos_unit     | String |    ""     | Unit in the DOS file | Extracted from the header if not set ``\\`` Currently supported: meV, eV, THz, Ry, Ha |
| nheader_dos  | Int64  |    -1     | Number of header lines in the DOS file | Auto-recognition if unset |
| nfooter_dos  | Int64  |    -1     | Number of footer lines in the DOS file | Auto-recognition if unset |

!!! details "Detailed description"
    >  **dos_file** :: STRING
    >
    > Path to the DOS file
    > Only required for vDOS calculations

    >  **spinDos** :: INTEGER | Default: 1
    >
    > Spin convention used in the DOS file
    >
    > | Value | Description |
    > | ---   | ----------- |
    > | 1     | Spin is not considered in the DOS |
    > | 2     | DOS contains spin degeneracy |

    >  **dos_unit** :: STRING    | Default: nothing
    >
    > Energy unit in the DOS file
    >
    > | Value | Description |
    > | ---   | ----------- |
    > | -1   | Auto-extraction from header |
    > | meV     | - |
    > | eV     | - |
    > | THz     | - |
    > | Ry     | - |
    > | Ha     | - |

    >  **nheader_dos** :: INTEGER | Default: nothing
    >
    > Number of header lines in the DOS file
    >
    > | Value | Description |
    > | ---   | ----------- |
    > | -1   | Auto-recognition |

    >  **nfooter_dos** :: INTEGER  | Default: nothing
    >
    > Number of footer lines in the DOS file
    >
    > | Value | Description |
    > | ---   | ----------- |
    > | -1   | Auto-recognition |


### Weep
The *Weep_file* is required for ``W`` calculations.
By default, IsoME assumes that the third column contains the ``W(\varepsilon,\varepsilon')`` values and that the first and second columns contain the corresponding energy-grid coordinates.
The columns can be changed via *Weep_col* and *Wen_col*.
If the row and column identifiers are consecutive indices rather than energies, provide an additional *Wen_file* containing the ``W`` energy-grid points.
#### Summary formatting:
- **header:** Non-numeric rows at the beginning of the document. If the header contains the unit (meV, eV, THz, Ry, Ha), it is extracted automatically; otherwise, set the unit via *Weep_unit*. If the header contains only one numeric value, it is interpreted as the Fermi energy. If there are several numerical values, IsoME checks for a keyword (`ef`, `efermi`, ...) indicating the Fermi energy. If extraction fails, adapt the header or set the Fermi energy via *efW*.
- **footer:** Non-numeric rows at the end of the document.
- **columns:** energy-grid points are assumed to be in the first and second columns, and the first column is used by default (*Wen_col* = 1). The ``W`` values are assumed to be in the third column by default (*Weep_col* = 3).


If the *Weep_file* does not contain the energies, an additional *Wen_file* can be specified.
- **header:** Non-numeric rows at the beginning of the document. It is assumed that the header contains the unit. If not, the unit has to be specified via *Wen_unit*.
- **footer:** Non-numeric rows at the end of the document.
- **first column:** energy grid points ``\epsilon`` of ``W(\epsilon, \epsilon')``. The column containing the energies can be changed via *Wen_col*.


|     Name     |  Type  |  Default  |          Description          |                Comment              |
|--------------|--------|:---------:|---------------------------|-----------------------------------------|
| Weep_unit    | String |     ""    | Unit of ``W``             | Extracted from the header of `Weep_file` if not set ``\\`` Currently supported: meV, eV, THz, Ry, Ha |
| Wen_unit     | String |     ""    | Unit of the ``W`` energies    | Extracted from the header of `Wen_file` if not set ``\\`` Currently supported: meV, eV, THz, Ry, Ha |
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

    >  **Weep_unit** :: STRING | Default: nothing
    >
    > Energy unit in the ``W`` file.
    >
    > | Value | Description |
    > | ---   | ----------- |
    > | nothing | Auto-extraction from header |
    > | meV     | - |
    > | eV      | - |
    > | THz     | - |
    > | Ry      | - |
    > | Ha      | - |

    >  **Wen_unit** :: STRING | Default: nothing
    >
    > Energy unit in the optional ``W`` energy-grid file.
    >
    > | Value | Description |
    > | ---   | ----------- |
    > | nothing | Auto-extraction from header |
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



# Version
Julia 1.10 or higher is required
