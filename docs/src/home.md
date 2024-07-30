# isoME
The package isoME.jl offers a quick way to solve the isotropic Eliashberg equations.  
The superconducting critical temperature (Tc) as well as the components of the self-energy $\Delta, Z, \chi$ are returned. (cite our paper??)  
In the simplest case only a file containing the Eliashberg spectral function $\alpha^2 F$ is required as input.
For more advanced calculations a file containing the DOS and the screened Coulomb interaction W, respectively are needed.
The package is able to solve the isotropic eliashberg equations within one of the following approximations:
- constant Dos + Anderson pseudopotential mu*
- full bandwidth + Anderson pseudopotential mu*
- constant Dos + screened Coulomb interaction W
- full bandwidth + screened Coulomb interaction W

Furthermore, results based on the Allen-Dynes and Allen-Dynes-McMillan formulas are provided. 

## Installation
 -------------------- to do ---------------------  
isoME.jl is (not yet) a registered Julia package and can thus be installed via

using Pkg
Pkg.add("isoME")

 ## Quick start
- a2f_file = "path to alpha2F"
- dos_file = "path to dos" [optional]  
    - spinDos = {1,2}  
    - colFermi_dos = column of fermi energy in header, starting with 0 for last column  
    - cDOS_flag = {0,1} 
- Weep_file = "path to W" [optional] 
    - dosW_file =  "path to W dos"
    - include_Weep = {0,1}     
- temps = Temperature range [optional]


 ## Outline ??

