# FAQ
## Q1: The conversion of the input files is wrong?
Check if the header of the input file contains an unexpected unit. Currently supported are **eV, meV, THz, Ry and Ha**. The units are case sensitive! 
If it is still not working enter the unit manually via the input parameters (a2f_unit, dos_unit, Weep_unit, Wen_unit)

## Q2: Why does the order parameter start to oscillate between positive and negative values?
The eliashberg equations are unstable when the coulomb part of the order parameter exceeds the phonon part at any iteration. Double check if the *muc_ME* is set correctly. Try setting the damping of the coulomb part to a higher value via *nItFullCoul* or use a different mixing factor. It is also possible that the material is simply not a superconductor at the given temperature.

## Q3: What does -1 and "" indicate in the input structure?
-1 or "" indicate that either a default value is used or that it is an optional input. If a mandatory input parameters is set to -1 or "" it will be overwritten during the run.

E.g. If the fermi-energy (ef) is unspecified its value will be extracted from the header of the dos-file and the -1 will be overwritten with the acutal value.

## Q4: Why is the code unable to read my input files?
Per default, a certain structure of the input files is assumed:
### a2F-file: 
- 1st column: energies
- 2nd - nth column: a2F-values for different smearings. Use *ind_smear* to select a specific column.
### dos-file:
- 1st column: energies
- 2nd column: dos
### Weep-file
- 1st column: energies (change via *Wen_col*)
- 3rd column: Weep data (change via *Weep_col*)
### Wen-file
Only required if the Weep file does not contain the energy grid points
- 1st column: energies (change via *Wen_col*)