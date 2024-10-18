# FAQ
## Q1: The conversion of the input files is wrong?
Check if the header of the input file contains an unexpected unit. Currently supported are **eV, meV, THz and Ry**. The units are case sensitive! 
If it is still not working enter the unit manually via the input parameters (a2f_unit, dos_unit, Weep_unit, Wen_unit)

## Q2: Why does the order parameter start to oscillate between positive and negative values?
The eliashberg equations are unstable when the coulomb part is greater than the phonon part at any iteration. Try setting the damping of the coulomb part to a higher value via nItFullCoul or using a different mixing factor. It is also possible that your muc_ME is too large or that the material is simply not a superconductor at the given temperature.

## Q3: What does -1 indicate in the input structure?
-1 indicates that either a default value is used or that it is an optional input. The value of all mandatory input parameters will be overwritten during the run.

E.g. If the fermi-energy (ef) is unspecified its value will be extracted from the header of the dos-file the -1 will be overwritten with the acutal value.

## Q4: Why is the code unable to read my input files?
Per default, a certain structure of the input files is assumed:
### a2F-file: 
- 1st column: energies
- 2nd - nth column: a2F-values for different smearings
### dos-file:
- 1st column: energies
- 2nd column: dos
### Weep-file
- 3rd column: Weep data (change via Weep_col)
### Wen-file
- 1st column: energies