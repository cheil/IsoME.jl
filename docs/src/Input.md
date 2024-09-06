# Input Parameters
The input parameters are handed over collectively via a structure (arguments). We recommend to always initialzie a new instance of the struct when running the Eliashberg solver as some of the parameters may be overwritten during a run.
All energies internally are assumed to be in meV.  
If the input files differ from that, the units are extracted from the header and converted automatically. If this doesn't work for some reason, the user has to specify the units manually via the corresponding input parameters.  
A comprehensive description of all the input paramters can be found below.

## General
General input parameters are:

| Name    |      Type      | Description | Comment  |
|---------|----------------|-------------|----------|
| temps   | Vector{Number} | Temperatures| Only needed if TcSearchMode_flag = 0|
| mu      | Float64        | Parameter for mu*, if specified muc_AD & muc_ME are calculated based on mu | |
| ef      | Float64        | Fermi-energy, in case mu is specified but no dos-file |  |
| muc_AD     | Float64        | Morel-Anderson Pseudopotential Allen-Dynes [cite] | Default = 0.14|
| muc_ME  | Float64        | Morel-Anderson Pseudopotential Migdal-Eliashberg [cite] | Default = muc / (1 + muc*log(omega_ph/omega_c)) , ref paper? |
| omega_c | Float64        | Matsubara cutoff in meV | |
| mixing_beta | Number     | Linear mixing factor | Default: Iteration dependent |
| nItFullCoul | Number     | First iteration in which full coulomb interaction is used , dampens oscillations | Default = 5 |
| cDOS_flag | Int64 | 0: constant dos , 1: variable dos  | Default = 1 |
|  include_Weep | Int64 | 0: Morel-Anderson Pseudopotential , 1: W | Default = 0 |
|TcSearchMode_flag | Int64 | 0: Manual mode (Loop over all values in temps) , 1: Automatic Tc search | Default = 1 |
| mu_flag | Int64 | Only if cDOS_flag = 0 , 0: no mu-update , 1: mu-update (recommended) | Default = 1 |
| outdir | String | path to the output directory | Default = pwd()*"/output/" |
| flag_log | Int64 | 0: no log-file , 1: create log-file | Default = 1 |
| flag_figure | Int64 | 0: no figures are plotted , 1: plot gap and a2f-values | Default = 1 |
| flag_outfile | Int64 | Create an output file | Default = 1 |
| material | String | Name of compound | Default = "Material" |    


Depending on the mode chosen, the following input files may be needed.

## $\alpha^2F$
A file containing the Eliashberg spectral function $\alpha^2F(\omega)$ is required for all calculations.  
It must obey the following structure:
- **header:** Optional. If the unit (meV, THz, Ry) is available it will be extraced automatically, otherwise set the unit via dos_unit. 
- **first column:** energies
- **second column onwards:** $\alpha^2F$ values for different smearings. Per default the smearing in the middle is used.
- **footer:** Optional

The number of header/footer lines and smearing values should be recognized automatically.

### Input
| Name        | Type | Description                           | Comment                       |
|-------------|----------|-----------------------------|-------------------------------|
| a2f_file    | String | path to $\alpha^2F$-file              | required for all calculations |
| ind_smear   | Int64 | index of smearing that should be used | optional                      |
| a2f_unit    | String| unit in which the data is given       | optional                      |
| nsmear      | Number | number of smearings in the a2f-file   | optional                      |
| nheader_a2f | Number | number of header lines in a2f_file    | optional                      |
| nheader_a2f | Number| number of footer lines in a2f_file    | optional  |

<details>
<summary>Detailed description</summary>
<br>

>  **a2f_file** :: STRING 
> 
> Path to the $\alpha^2F$-file.  
> Required for all calculations.  
> The first column must contain the energies, the second column onwards the $\alpha^2F$ values for different smearings.  
> 

> **ind_smear** :: INTEGER | Default: 1
>
> Index of the smearing which should be used.
> 

>  **a2f_unit** :: STRING    | Default: nothing
>
> Energy unit in the $\alpha^2F$-file
>
> | Value | Description |
> | ---   | ----------- |
> | nothing   | Auto-extraction from header |
> | meV     | - |
> | THz     | - |
> | Ry     | - |

>  **nheader_a2f** :: INTEGER | Default: nothing  
> 
> Number of header lines in the $\alpha^2F$-file
>
> | Value | Description |
> | ---   | ----------- |
> | nothing   | Auto-recognition |

>  **nfooter_a2f** :: INTEGER  | Default: nothing  
>
> Number of footer lines in the $\alpha^2F$-file
>
> | Value | Description |
> | ---   | ----------- |
> | nothing   | Auto-recognition |
</details>                    

## $N(\epsilon)$
For a full bandwidth calculation a file containing the dos is required.  
It must obey the following structure:
- **header:** It is assumed that the fermi energy is the second to last entry of the first line as in a QE Dos-file. If this is not the case, specify the position of the fermi-energy via colFermi_dos starting from the right with 0 for the last column. Furthermore, it is assumed that the header contains the unit. If not, the unit has to be specified via dos_unit
- **first column:** energies
- **second column:** dos values
- **footer:** Optional 

| Name        | Description                     | Comment                                                |
|-------------|---------------------------------|--------------------------------------------------------|
| **dos_file**    | path to the dos-file            | Only required for FBW calculations 
| **colFermi_dos** | column of the fermi energy starting from the right | 0 = last column, 1 = second to last column, ... |                    |
| **spinDos**     | double count of dos due to spin | 1 = no double count, 2 = double count                  |
| **dos_unit**    | unit of the energies            | Optional, auto-extraction from header possible         |
| **nheader_dos** | number of header lines          | Optional, auto recognition if the header contains text |
| **nfooter_dos** | number of footer lines          | Optional, auto recognition if the footer contains text |
|             |                                 |                                                        |
<details>
<summary>Detailed description</summary>
<br>

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

## Weep
For a calculation using the full screened Coulomb interaction two files containing $W(\epsilon,\epsilon')$ and the energy grid $\epsilon$, respectively are needed in addition to the dos-file.  
The $W$-file must obey the following structure:
- **header:** It is assumed that the header contains the unit. If not, the unit has to be specified via Weep_unit
- **column:** Specify column the W data is in via Weep_col (Default=3)
- **footer:** Optional 

The energy-file must obey the following structure:
- **header:** It is assumed that the header contains the unit. If not, the unit has to be specified via Wen_unit
- **first column:** energy grid points $\epsilon$ of $W(\epsilon, \epsilon')$
- **footer:** Optional

### Input
| Name        | Type | Description                           | Comment                       |
|-------------|----------|-----------------------------|-------------------------------|
| Weep_file | String | Path to W-file | |
| Wen_file| String | Path to file containing the W energy grid points | |
| Weep_unit | String | unit in which the W data is given  |optional |                 
| Weep_col | Number | column the W data is in | Default = 3 |
| nheader_Weep | Number | number of header lines in the W-file | Default = Autoextraction |
| nfooter_Weep | Number | number of footer lines in the W-file | Default = Autoextraction |
| Wen_unit | String | unit in which the energies are given  |optional |              
| nheader_Wen | Number | number of header lines in the W Energy-file | Default = Autoextraction |
| nfooter_Wen | Number | number of footer lines in the W Energy-file | Default = Autoextraction |

<details>
<summary>Detailed description</summary>
<br>

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

# Examples

# Output

# Version
this package requires Julia 1.10 or higher (printf)