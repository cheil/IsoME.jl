# FAQ
## Q1: The conversion of the input files is wrong?
Check if the header of the input file contains an unexpected unit. Currently supported are **eV, meV, THz, Ry and Ha**. The units are case sensitive! 
If it is still not working enter the unit manually via the input parameters (a2f\_unit, dos\_unit, Weep\_unit, Wen\_unit)

## Q2: Why does the order parameter start to oscillate between positive and negative values?
The eliashberg equations are unstable when the coulomb part of the order parameter exceeds the phonon part at any iteration. Double check if the *muc_ME* is set correctly. Try setting the damping of the coulomb part to a higher value via *nItFullCoul* or use a different mixing factor. It is also possible that the material is simply not a superconductor at the given temperature.

## Q3: What does `NaN`, -1 and "" indicate in the input structure?
These sentinel values indicate that either a default value is used or that the input is optional: `NaN` for real-valued fields, `-1` for integer fields (line counts, column indices) and `""` for strings. If such an input is left at its sentinel, it will be determined during the run and the placeholder will be overwritten.

E.g. if the Fermi energy (*ef*) is unspecified (`NaN`), its value will be extracted from the header of the dos-file and the `NaN` will be replaced by the actual value.

## Q4: Why is the code unable to read my input files?
Per default, a certain structure of the input files is assumed:
#### a2F-file: 
- 1st column: energies
- 2nd - nth column: a2F-values for different smearings. Use *ind_smear* to select a specific column.
#### dos-file:
- 1st column: energies
- 2nd column: dos
#### Weep-file
- 1st column: energies (change via *Wen_col*)
- 3rd column: Weep data (change via *Weep_col*)
#### Wen-file
Only required if the Weep file does not contain the energy grid points
- 1st column: energies (change via *Wen_col*)

## Q5: My output directory only contains a `log.txt` — did the run fail?
It did not *crash*: a crash always leaves a `CRASH` file behind. A lone `log.txt` means the run
never reached its last step, either because it is still running — `Summary.dat` is only written at
the very end — or because it was stopped, by `Ctrl+C` or by something killing the process. The log
looks the same in all of these cases, so check whether the calculation is still running. See
[The output directory contains only a log.txt](@ref) on the [Troubleshooting](@ref) page.

## Q6: Which role do the cutoffs and grids play in the real-axis solver?
The real-axis solver is far more sensitive to its grids than the imaginary-axis solver: its integrands are complex, sharply peaked and their poles move during the iteration.
Most of the failures reported so far can be traced back to *omega_c*, *depsilon* or *domega*.
Each of these is discussed, together with the symptoms it produces, in the [Real-axis grids and the kernel](@ref) and [The μ-update](@ref) sections of the [Troubleshooting](@ref) page.
A detailed description of where each parameter enters is given on the [Real Axis Solver](@ref) page.