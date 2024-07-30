# Getting started
All energies internally are assumed to be in meV.  
If the input files differ from that, the units are extracted from the header and converted automatically. If this doesn't work for some reason, the user has to specify the units manually via the corresponding input parameters.  
A comprehensive description of all the input paramters can be found below.
## $\alpha^2F$
The Eliashberg spectral function $\alpha^2F(\omega)$ is required for all calculations.
### Input
| Name        | Description                           | Comment                       |
|-------------|---------------------------------------|-------------------------------|
| a2f_file    | path to $\alpha^2F$-file              | required for all calculations |
| ind_smear   | index of smearing that should be used | optional                      |
| a2f_unit    | unit in which the data is given       | optional                      |
| nsmear      | number of smearings in the a2f-file   | optional                      |
| nheader_a2f | number of header lines in a2f_file    | optional                      |
| nheader_a2f | number of footer lines in a2f_file    | optional  |

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
The density of states is required for the full bandwidth calculations.  
The first  column must contain the energies and the second column the corresponding dos values.
Furthermore, the file must have at least one header line containing at least the fermi energy.  
If the Fermi energy is not the last element of the first line, one has to specify the column it is in via colFermi_dos starting from the right with 0 for the last column.  
The unit is extracted from the header automatically if available. 
### Input
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
File containing the screened coulomb interactions.   
Furthermore, a file containg the dos on the energy grid of the W is needed.   
The other dos-file is optional in this case.
### Variables
- Weep_file: path to file containing the screened coulomb interaction
    * The values for the screened coulomb interaction must be in the third column


# Examples

# Output

# Version
this package requires Julia 1.10 or higher (printf)