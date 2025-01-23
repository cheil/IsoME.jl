# Input Parameters
The input parameters are handed over collectively as [compsite type](https://docs.julialang.org/en/v1/manual/types/#Composite-Types) with the name `arguments()`.  
We highly recommend to always initialzie a new instance of this struct when running the Eliashberg solver, because some of the parameters may be overwritten during a run.  
All energies internally are assumed to be in meV.  
If the input files differ from that, the units are extracted from the header and converted automatically. If this doesn't work for some reason, please double check the header of your input files or specify the units manually via the corresponding input parameters.  
A comprehensive description of all the input parameters can be found below.

## General
General input parameters are:

| Name    |      Type      |   Default   | Comment  | Description | 
|:--------|:---------------|:------------|:---------|:----------  |
| temps   | Vector{Number} |  [-1] | -1: Automatic Tc search mode ``\\`` else: Tc search at specified values|  Considered temperatures | 
| a2f_file    | String   |    ""    | Required for all calculations | Path to ``\alpha^2F``-file   
| ind_smear   | Int64    |    -1    | Per default the smearing in the middle column is used | Which smearing should be used | 
| a2f_unit    | String   |    ""    | Is extracted from the header of the dos-file if not set ``\\`` Currently supported: meV, eV, THz, Ry | unit of the a2f-file | 
| omega_c | Float64        | 7000 | Matsubara cutoff ``\omega_c`` in meV | |
| mixing_beta | Number     | Iteration dependent | Linear mixing factor |   |
| cDOS_flag | Int64 |  1   |0: variable dos <\br> 1: constant dos  | dos_file has to be specified |
| dos_file  |  String  |     ""    | Required cDOS_flag = 0| path to the dos-file |  
| ef        | Float64  |     -1    | Is extracted from the header of the dos-file if not set | Fermi-energy DOS in meV |
| colFermi_dos | Int64 |      1    |  0 = last column ``\\`` 1 = second to last column ``\\`` 2 = ... |       column the fermi energy is in starting from the right |
| spinDos   | Int64    |      2    |  1 = spin not considered ``\\`` 2 = spin considered  | Does the dos consider spin | 
| dos_unit   |  String |     ""    |  Is extracted from the header of the dos-file if not set ``\\`` Currently supported: meV, eV, THz, Ry | units in dos file | 
| mu      | Float64        |   -1  | Measure for the Coulomb strength (This is not the chemical potential!) ``\\`` If muc_AD & muc_ME are unspecified they are calculated based on mu | ``\mu=N(e_f)*W(e_f,e_f)`` ``\\``  cite paper with formula |
| muc_AD     | Float64     |  0.12 | Morel-Anderson Pseudopotential Allen-Dynes ``\mu^*_{AD}`` | cite paper with formula |
| muc_ME  | Float64        | ``{\mu^*_{AD}}/({1 + \mu^*_{AD} \ln(\omega_{ph}/\omega_c)})`` | Morel-Anderson Pseudopotential Migdal-Eliashberg ``\mu^*_{ME}`` | cite paper with formula |
| nItFullCoul | Number     |  5    | First iteration in which full coulomb interaction is used , dampens oscillations |  |
| mu_flag   | Int64 |  1   |  0: no mu-update ``\\`` 1: mu-update (recommended) | Update the chemical potential in vDOS calculations |
| include_Weep | Int64 | 0 |0: Morel-Anderson Pseudopotential ``\\`` 1: static Coulomb interaction W(e,e') | Weep_file and Wen_file have to be specified |
| Weep_file |  String  |     ""    |  Required if include_Weep = 1 |   Path to W-file  |
| Wen_file  |  String  |     ""    |  Required if include_Weep = 1 | Path to file containing the energy grid points of W |
| efW       | Float64  |     -1    | Is extracted from the header of the Wen-file if not set | Fermi-energy W energies in meV |
| Weep_unit |  String  |     ""    | Is extracted from the header of the Weep-file if not set ``\\`` Currently supported: meV, eV, Ry, Ha | Unit of the W | 
| Wen_unit  |  String  |     ""    | Is extracted from the header of the Wen-file if not set ``\\`` Currently supported: meV, eV, Ry, Ha  | Unit of the W energies |                 
| Weep_col  |   Int64  |     3     |                               | Column of the W-data in the Weep-file  | 
| outdir | String |  pwd() | Path to the output directory |  |
| flag_figure | Int64 |  1 | 0: no figures are plotted ``\\`` 1: plot gap and a2f-values |  |
| flag_writeSelfEnergy | Int64 | 0  | 0: Don't save the self-energy components 1: save the self-energy components at each temperature  | Specify if the self-energy components should be saved |
| material | String | "Material" | Title used in plots, summary, ... | Name of compound  |    
| encut  | Float64 |  -1 (full window is used) | in meV | Considered energy window around the fermi energy used in the integration of phic |
| shiftcut  | Float64 | 2000 meV  | in meV | Considered energy window around the fermi energy used in the integration of the shift. Alwas smaller than encut|


### ``\mu, \mu^*_{AD} \& \mu^*_{ME}``
``\mu`` measures the strength of the Coulomb interaction at the fermi-surface: ``\mu=W(\varepsilon_F,\varepsilon_F)``   
It is connected to the pseudopotentials via:   
``\mu*_{AD}=\frac{\mu}{1+\mu \text{ ln}\left(\frac{\varepsilon_{el}}{\hbar \omega_{ph}}\right)}``
where ``\omega_{ph}`` is a characteristic cutoff frequency for the phonon-induced interaction and ``\varepsilon_{el}`` is a characteristic electronic energy scale.   
A different phonon cutoff has to be used for the Allen-Dynes and Migdal-Eliashberg pseuodpotential.   
If the user does specify one of the ``\mu's`` all other will be calculated based on it, but the code will never overwrite an user input.  
Furthermore, if none of the ``\mu's`` is specified but a Weep-file exists, the ``\mu`` will be calculated from the ``W``.   
If the user specifies nothing, a default value of ``\mu^*_{AD}=0.12`` will be used. 
For the conversion the user can specify the typical electronic energy explicitly (input: typEl), otherwise the fermi energy (input: ef or efW)  will be used. For the characteristic phonon cutoff the Matsubara cutoff  or the maximum given phonon frequency will be used in case of ME or AD, respectively.

## Expert user input parameters
It should not be necessary to change any of the following input parameters during a normal execution of the code.
They are only relevant if you encounter any problems.
|     Name     |  Type  |  Default  |          Comment          |                Description              | 
|--------------|--------|:---------:|---------------------------|-----------------------------------------|
| nsmear       | Int64  |    -1     | Auto-recognition if unset | number of smearings in the a2f-file     |
| nheader_a2f  | Int64  |    -1     | Auto-recognition if unset | number of header lines in a2f_file      | 
| nheader_a2f  | Int64  |    -1     | Auto-recognition if unset | number of footer lines in a2f_fi        |
| nheader_dos  | Int64  |    -1     | Auto-recognition if unset | number of header lines in dos-file      |
| nfooter_dos  | Int64  |    -1     | Auto-recognition if unset | number of footer lines in dos-file      | 
| nheader_Weep | Int64  |    -1     | Auto-recognition if unset | number of header lines in the Weep-file |
| nfooter_Weep | Int64  |    -1     | Auto-recognition if unset | number of footer lines in the Weep-file |          
| nheader_Wen  | Int64  |    -1     | Auto-recognition if unset | number of header lines in the Wen-file  |
| nfooter_Wen  | Int64  |    -1     | Auto-recognition if unset | number of footer lines in the Wen-file  |


## Input files
Depending on the mode chosen, the following input files may be needed.

### ``\alpha^2F``
A file containing the Eliashberg spectral function ``\alpha^2F(\omega)`` is required for all calculations.  
It must obey the following structure:
- **header:** Optional. If the unit (meV, THz, Ry) is available it will be extraced automatically, otherwise set the unit via dos_unit. 
- **first column:** energies
- **second column onwards:** ``\alpha^2F`` values for different smearings. Per default the smearing in the middle is used.
- **footer:** Optional

The number of header/footer lines and smearing values should be recognized automatically.

<details>
<summary>Detailed description</summary>
``\\``

>  **a2f_file** :: STRING 
> 
> Path to the ``\alpha^2F``-file.  
> Required for all calculations.  
> The first column must contain the energies, the second column onwards the ``\alpha^2F`` values for different smearings.  
> 

> **ind_smear** :: INTEGER | Default: 1
>
> Index of the smearing which should be used.
> 

>  **a2f_unit** :: STRING    | Default: nothing
>
> Energy unit in the ``\alpha^2F``-file
>
> | Value | Description |
> | ---   | ----------- |
> | nothing   | Auto-extraction from header |
> | meV     | - |
> | THz     | - |
> | Ry     | - |

>  **nheader_a2f** :: INTEGER | Default: nothing  
> 
> Number of header lines in the ``\alpha^2F``-file
>
> | Value | Description |
> | ---   | ----------- |
> | nothing   | Auto-recognition |

>  **nfooter_a2f** :: INTEGER  | Default: nothing  
>
> Number of footer lines in the ``\alpha^2F``-file
>
> | Value | Description |
> | ---   | ----------- |
> | nothing   | Auto-recognition |
</details>                    

### ``N(\epsilon)``
For a full bandwidth calculation a file containing the dos is required.  
It must obey the following structure:
- **header:** It is assumed that the fermi energy is the second to last entry of the first line as in a QE Dos-file. If this is not the case, specify the position of the fermi-energy via colFermi_dos starting from the right with 0 for the last column. Furthermore, it is assumed that the header contains the unit. If not, the unit has to be specified via dos_unit
- **first column:** energies
- **second column:** dos values
- **footer:** Optional 

<details>
<summary>Detailed description</summary>
``\\``

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
</details>

### Weep
For a calculation using the full screened Coulomb interaction two files containing ``W(\epsilon,\epsilon')`` and the energy grid ``\epsilon``, respectively are needed in addition to the dos-file.  
The ``W``-file must obey the following structure:
- **header:** It is assumed that the header contains the unit. If not, the unit has to be specified via Weep_unit
- **column:** Specify column the W data is in via Weep_col (Default=3)
- **footer:** Optional 

The energy-file must obey the following structure:
- **header:** It is assumed that the header contains the unit. If not, the unit has to be specified via Wen_unit
- **first column:** energy grid points ``\epsilon`` of ``W(\epsilon, \epsilon')``
- **footer:** Optional


<details>
<summary>Detailed description</summary>
``\\``

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
</details>
</br>



# Version
Julia 1.10 or higher is required (printf)