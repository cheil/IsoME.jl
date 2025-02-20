# Input
The input parameters are handed over collectively as [compsite type](https://docs.julialang.org/en/v1/manual/types/#Composite-Types) with the name `arguments()`.  
We highly recommend to always initialzie a new instance of this struct when running the Eliashberg solver, as some of the parameters may be overwritten during a run. Such parameters are marked by either `-1` or `""` for numeric and string variables, respectively.  

All energies internally are assumed to be in meV. If the input files differ from that, the units are extracted from the header and converted automatically. If this doesn't work for some reason, please double check the header of your input files or specify the units manually via the corresponding input parameters.  
A comprehensive description of all the input parameters can be found below.

## General
An overview of the most relevant inputs is given below. Keep in mind that all `-1`'s or `""`'s of required variables will be overwritten during the execution.  For more details please refer to the dedicated sections below.

 Name    |      Type      |   Default   | Description | Comment  |
|:--------|:---------------|:------------|:------------|:---------|
| temps   | Vector{Number} |  [-1] | Considered temperatures | -1: Tc search mode ``\\`` else: Tc search at specified values |
| a2f_file    | String   |    ""    | Path to ``\alpha^2F``-file | Required for all calculations |
| ind_smear   | Int64    |    -1    | Used smearing column  | Per default the smearing in the middle column is used |
 mu      | Float64        |   -1  | ``\mu=N(e_f)*W(e_f,e_f)`` | Measure for the Coulomb strength |
| muc_AD     | Float64     |  -1 | Morel-Anderson Pseudopotential Allen-Dynes ``\mu^*_{AD}`` | - |
| muc_ME  | Float64        | -1 | Morel-Anderson Pseudopotential Migdal-Eliashberg ``\mu^*_{ME}`` | - |
| typEl | Float64|  -1 | typical electronic energy |  used to calculate ``\mu^*`` from ``\mu`` |
| omega_c | Float64        | 7000 | Matsubara cutoff ``\omega_c`` in meV | |
| encut  | Float64 |  5000  | Cutoff for integration | in meV |
| shiftcut  | Float64 | 2000   | Cutoff for integration of shift. Always smaller than encut | in meV |
| mixing_beta | Number     | Iteration ``\\`` dependent | Linear mixing factor | |
| cDOS_flag | Int64 |  1   | dos_file has to be specified | 0: variable dos ``\\`` 1: constant dos |
| dos_file  |  String  |     ""    | Path to the dos-file | Required if cDOS_flag = 0 |
| ef        | Float64  |     -1    | Fermi-energy DOS in meV | Is extracted from the header of the dos-file if not set |
| mu_flag    | Int64 | 1   |  Update the chemical potential in vDOS calculations? | 0: no ``//`` 1: yes (recommended) |
| nItFullCoul | Number     |  10    | Dampens Coulomb contribution unitl *nItFullCoul*th iteration | |
| include_Weep | Int64 | 0 | *Weep_file* and *Wen_file* have to be specified | 0: Morel-Anderson Pseudopotential ``\\`` 1: static Coulomb interaction ``W(\varepsilon,\varepsilon')`` |
| Weep_file |  String  |     ""    | Path to W-file | Required if *include_Weep* = 1 |
| Wen_file  |  String  |     ""    | Path to file containing the energy grid points of W | Only required if the energies are not contained in the *Weep_file*|
| efW       | Float64  |     -1    | Fermi-energy W in meV | Is extracted from the header of the Wen-file if not set |
| conv_thr | ``10^{-4}`` | convergence threshold | - |
| minGap   | Float64 |0.1  | termination criterion | ``\Delta(0)<`` *minGap* |
| N_it | Int64 | 5000 | maximum number of iterations | - |
| sparseSamplingTemp | Float64| 2 | maximum temperature for sparse sampling | - |
| itpBounds | Vector{Float64} | [100,500] | interpolation intervals | - |
| itpStepSize | Vector{Int64} | [1,5,50] | interpolation step size | - |
| outdir | String |  pwd() | Path to the output directory | |
| flag_figure | Int64 |  1 | Should the gap and the ``\alpha^2F``-values be plotted? | 0: No ``\\`` 1: Yes |
| flag_writeSelfEnergy | Int64 | 0  | Should the self-energy components be saved? | 0: No ``\\`` 1: Yes |
| material | String | "Material" | Name of compound | Title used in plots, summary, ... |




### Pseudopotentials ``\mu,~\mu^*_{AD}~\&~\mu^*_{ME}``
``\mu`` measures the strength of the Coulomb interaction at the fermi-surface: ``\mu=W(\varepsilon_F,\varepsilon_F)``. 
It is connected to the pseudopotentials via:   
```math
\mu^*_{AD}=\frac{\mu}{1+\mu \text{ ln}\left(\frac{\varepsilon_{el}}{\hbar \omega_{ph}}\right)}
```
where ``\omega_{ph}`` is a characteristic cutoff frequency for the phonon-induced interaction and ``\varepsilon_{el}`` is a characteristic electronic energy scale.   
The typical electronic energy can be explicitly sepcified through *typEl*, otherwise the fermi energy (*ef* or *efW*)  will be used. For the characteristic phonon cutoff the Matsubara cutoff and the maximum given phonon frequency are used in case of ME or AD, respectively.

Per default, ``\mu`` and ``\mu^*``'s are set to -1, which indicates that *muc_AD* = 0.12 is used. This also fixes *muc_ME* through 
```math
\mu^*_{ME}={\mu^*_{AD}}/({1 + \mu^*_{AD} \ln(\omega_{ph}/\omega_c)})~.
```
If at least one of the ``\mu's`` is set, all unspecified ``\mu^*``'s, will be calculated based on it. However, no user input will be overwritten.
Furthermore, if none of the ``\mu``'s is specified but a Weep-file is given, the ``\mu`` will be calculated from the ``W``.   


## Read-In 
IsoME automatically recognizes the formatting of ``\textsc{QE/EPW/BerkeleyGW}`` files. 
Compatibility with other DFT/DFPT/GW packages is currently under development. However, the read-in function has been designed to rely on as little formatting as possible and will often work for different formattings as well. If the auto-recognition fails, either adapt the formatting of the input files or manually set it through the dedicated flags. For details please refer to the dedicated section for each input file.
Possible sources of errors are the number of header/footer lines, the Fermi energy, or the units. 

In its auto-recognition mode, IsoME interprets non-numerical rows at the beginning (end) of the file as the header (footer). From the header, the unit (Currently supported units are meV, eV, THz, Ry, Ha) and, in case of a dos/W-file, the Fermi energy are extracted. 


### ``\alpha^2F``
The Eliashberg spectral function ``\alpha^2F(\omega)`` is required for all calculations.  
The ``\alpha^2F(\omega)`` -file can contain an arbitrary amount of columns, but the first column must contain the energies and the remaining columns are interpreted as ``\alpha^2F(\omega)`` -values for different smearings. If the user does not specify the smearing column via *ind_smear*, the column in the middle will be used.


#### Summary formatting:
- **header:** Non-numeric rows at the beginning of the document. If the header contains the unit (meV, eVm THz, Ry, Ha), it will be extraced automatically, otherwise set the unit via *dos_unit*.  
- **footer:** Non-Numeric rows at the end of the document.
- **first column:** energies
- **second column onwards:** ``\alpha^2F`` values for different smearings. Per default the smearing in the middle is used.

The number of header/footer lines and smearing values should be recognized automatically. If this is not the caseset it through the dedicated input parameters:

|     Name     |  Type  |  Default  |     Description       |     Comment         | 
|--------------|--------|:---------:|-------------------------------------|---------------------------|
| a2f_unit     | String |    ""     | Unit of the a2f-file                | Is extracted from the header of the dos-file if not set. Currently supported: meV, eV, THz, Ry, Ha |
| nsmear       | Int64  |    -1     | number of smearings in the a2f-file | Auto-recognition if unset | 
| nheader_a2f  | Int64  |    -1     | number of header lines in a2f_file  | Auto-recognition if unset | 
| nheader_a2f  | Int64  |    -1     | number of footer lines in a2f_fi    | Auto-recognition if unset | 

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
For vDOS calculations a dos-file is required.  
The first and second column of the dos-file are interpreted as the energies and dos values, respectively. All other columns are ignored. The dos-values are divided by two, to remove the double counting due to spin. If this is not desired set the *spinDos* flag to 1. 

#### Summary formatting:
- **header:** Non-numeric rows at the beginning of the document. If the header contains the unit (meV, eVm THz, Ry, Ha), it will be extraced automatically, otherwise set the unit via *dos_unit*. If the header contains only one numeric value this is interpreted as the fermi energy. If there are several numerical values, IsoME checks if a keyword (ef, efermi,...) indicating the fermi energy exists. If the extraction of the fermi energy fails either adapt your header or set the fermi energy via *ef*.
- **footer:** Non-Numeric rows at the end of the document.
- **first column:** energies
- **second column:** dos values

|     Name     |  Type  |  Default  |          Description          |                Comment              | 
|--------------|--------|:---------:|---------------------------|-----------------------------------------|
| spinDos      | Int64  |     2     | Does the dos consider spin | 1 = spin not considered ``\\`` 2 = spin considered |
| dos_unit     | String |    ""     | Units in dos file         | Is extracted from the header of the dos-file if not set ``\\`` Currently supported: meV, eV, THz, Ry, Ha |
| nheader_dos  | Int64  |    -1     | number of header lines in dos-file      |  Auto-recognition if unset |
| nfooter_dos  | Int64  |    -1     |  number of footer lines in dos-file      | Auto-recognition if unset |

!!! details "Detailed description"
    >  **dos_file** :: STRING 
    > 
    > Path to the dos-file  
    > Only required for vDOS calculations

    >  **spinDos** :: INTEGER | Default: 1
    >
    > Double count of dos due to spin
    >
    > | Value | Description |
    > | ---   | ----------- |
    > | 1     | Spin is not considered in the dos |
    > | 2     | Double count of dos due to spin |

    >  **dos_unit** :: STRING    | Default: nothing
    >
    > Energy unit in the dos file
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
    > Number of header lines in the dos file
    >
    > | Value | Description |
    > | ---   | ----------- |
    > | -1   | Auto-recognition |

    >  **nfooter_dos** :: INTEGER  | Default: nothing  
    >
    > Number of footer lines in the dos file
    >
    > | Value | Description |
    > | ---   | ----------- |
    > | -1   | Auto-recognition |


### Weep
The *Weep_file* is required for W calculations.
It is assumed that the third row contains the $W(\varepsilon,\varepsilon')$ values and the first and second columns the row and column numbering as energy values. The respective columns can be change via *Weep_col* and *Wen_col*.
 If the row/column numbering is given as consecutive numbering, an additional input file, called *Wen_file*, containing the $W$ energy grid points, can be handed over. 
#### Summary formatting:
- **header:** Non-numeric rows at the beginning of the document. If the header contains the unit (meV, eVm THz, Ry, Ha), it will be extraced automatically, otherwise set the unit via *Weep_unit*. If the header contains only one numeric value, it is interpreted as the fermi energy. If there are several numerical values, IsoME checks if a keyword (ef, efermi,...) indicating the fermi energy exists. If the extraction of the fermi energy fails either adapt your header or set the fermi energy via *efW*.
- **footer:** Non-Numeric rows at the end of the document. 
- **column:** energy-grid points are assumed to be in the first and second column and the first column is used per default (*Wen_col*=1). The W-values are assumed to be in the third column per default (*Weep_col*=3).


If the *Weep_file* does not contain the energies an additional *Wen_file* can be specified.
- **header:** Non-numeric rows at the beginning of the document. It is assumed that the header contains the unit. If not, the unit has to be specified via *Wen_unit*.
- **footer:** Non-Numeric rows at the end of the document. 
- **first column:** energy grid points ``\epsilon`` of ``W(\epsilon, \epsilon')``. The column containing the energies can be changed via *Wen_col*.


|     Name     |  Type  |  Default  |          Description          |                Comment              | 
|--------------|--------|:---------:|---------------------------|-----------------------------------------|
| Weep_unit    | String |     ""    | Unit of the W             | Is extracted from the header of the Weep-file if not set ``\\`` Currently supported: meV, eV, THz, Ry, Ha |
| Wen_unit     | String |     ""    | Unit of the W energies    | Is extracted from the header of the Wen-file if not set ``\\`` Currently supported: meV, eV, THz, Ry, Ha  |
| Weep_col     | Int64  |     3     | Column of the W-data in the Weep-file | |
| nheader_Weep | Int64  |    -1     | number of header lines in the Weep-file |Auto-recognition if unset | 
| nfooter_Weep | Int64  |    -1     | number of footer lines in the Weep-file |  Auto-recognition if unset |         
| nheader_Wen  | Int64  |    -1     |  number of header lines in the Wen-file  | Auto-recognition if unset | 
| nfooter_Wen  | Int64  |    -1     | number of footer lines in the Wen-file  | Auto-recognition if unset | 

!!! details "Detailed description"
    >  **dos_file** :: STRING 
    > 
    > Path to the dos-file  
    > Only required for FBW calculations
    
    > **colFermi_dos** :: INTEGER | Default: 0
    >
    > The fermi-energy has to be part of the header in the dos file.  
    > Specify the column it is in where per convetion 0 denotes the last column, 1 the second to last column and so on.
    > 

    >  **spinDos** :: INTEGER | Default: 1
    >
    > Specify if dos is multiplied by 2 due to spin. 
    >
    > | Value | Description |
    > | ---   | ----------- |
    > | 1     | Spin is not considered in the dos |
    > | 2     | Double count of dos due to spin |

    >  **dos_unit** :: STRING    | Default: nothing
    >
    > Energy unit in the dos file
    >
    > | Value | Description |
    > | ---   | ----------- |
    > | nothing   | Auto-extraction from header |
    > | meV     | - |
    > | eV     | - |

    >  **nheader_dos** :: INTEGER | Default: nothing  
    > 
    > Number of header lines in the dos file
    >
    > | Value | Description |
    > | ---   | ----------- |
    > | nothing   | Auto-recognition |

    >  **nfooter_dos** :: INTEGER  | Default: nothing  
    >
    > Number of footer lines in the dos file
    >
    > | Value | Description |
    > | ---   | ----------- |
    > | nothing   | Auto-recognition |



# Version
Julia 1.10 or higher is required 