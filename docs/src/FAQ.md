# FAQ
## Q1: Why are my input files converted to the wrong units?
Check if the header of the input file contains an unexpected unit. Currently supported are **meV, eV, THz, Ry and Ha**. The units are case sensitive!
If no unit can be read from the header at all, the run stops with `Could not determine the unit of the …-file`, which names the input to set.
Either way, the unit can be given manually via `a2f_unit`, `dos_unit`, `Weep_unit` or `Wen_unit`.

## Q2: Why does the order parameter start to oscillate between positive and negative values?
The Eliashberg equations are unstable when the Coulomb part of the order parameter exceeds the phonon part at any iteration. Double check if `muc_ME` is set correctly. Try setting the damping of the Coulomb part to a higher value via `nItFullCoul`, or use a different mixing factor via `mixing_beta`. It is also possible that the material is simply not a superconductor at the given temperature.

## Q3: What does `NaN`, -1 and "" indicate in the input structure?
These sentinel values indicate that either a default value is used or that the input is optional: `NaN` for real-valued fields, `-1` for integer fields (line counts, column indices) and `""` for strings. If such an input is left at its sentinel, it will be determined during the run and the placeholder will be overwritten.

E.g. if the Fermi energy (`ef`) is unspecified (`NaN`), its value will be extracted from the header of the dos-file and the `NaN` will be replaced by the actual value.

Note that `-1` is a sentinel for integer fields only. Before IsoME 2.0 the real-valued fields used it as well, so an explicit `-1` on `mu`, `muc_AD`, `muc_ME`, `ef`, `efW`, `typEl` or `mixing_beta` is now rejected with a message pointing at `NaN`.
The same check rejects negative values for `mu`, `muc_AD`, `muc_ME`, `mixing_beta`, `conv_thr`, `minGap`, `N_it` and `min_it`.
See [Errors from `arguments()` itself](@ref) on the [Troubleshooting](@ref) page for the messages `arguments()` can produce.

## Q4: Why is the code unable to read my input files?
By default, a certain structure of the input files is assumed:
#### a2F-file:
- 1st column: energies
- 2nd - nth column: a2F-values for different smearings. Use `ind_smear` to select a specific column.
#### dos-file:
- 1st column: energies
- 2nd column: dos
#### Weep-file:
- 1st column: energies (change via `Wen_col`)
- 3rd column: Weep data (change via `Weep_col`)
#### Wen-file:
Only required if the Weep file does not contain the energy grid points.
- 1st column: energies (change via `Wen_col`)

## Q5: Why is my output not in the directory I gave?
An existing directory is never written into, so that a rerun cannot overwrite earlier results.
If `outdir` already exists, a run counter is appended instead: `outdir`, `outdir_1`, `outdir_2`, and so on.
With the default `outdir` this happens on every run after the first, and the results land in `IsoME/`, `IsoME_1/`, `IsoME_2/`, … inside the current working directory.

## Q6: My output directory only contains a `log.txt` — did the run fail?
It did not *crash*: a crash always leaves a `CRASH` file behind. A lone `log.txt` means the run
never reached its last step, either because it is still running — `Summary.dat` is only written at
the very end — or because it was stopped by the user (e.g. `Ctrl+C`). The log
looks the same in all of these cases, so check whether the calculation is still running. See
[The output directory contains only a log.txt](@ref) on the [Troubleshooting](@ref) page.

## Q7: Which role do the cutoffs and grids play in the real-axis solver?
The real-axis solver is far more sensitive to its grids than the imaginary-axis solver: its integrands are complex, sharply peaked and their poles move during the iteration.
Most of the failures reported so far can be traced back to `omega_c` or `encut`. Also `depsilon`, `domega` and `n_cheb` can play a role, although the default parameters have been sufficient to achieve convergence in our benchmarks.
Each of these is discussed, together with the symptoms it produces, in the [Real-axis grids and the kernel](@ref) and [The μ-update](@ref) sections of the [Troubleshooting](@ref) page.
A detailed description of where each parameter enters is given on the [Real Axis Solver](@ref) page.
