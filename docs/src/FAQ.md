# FAQ
## Q1: The conversion of the input files is wrong?
Check if the header of the input file contains an unexpected unit. Currently supported are **eV, meV, THz, Ry and Ha**. The units are case sensitive! 
If it is still not working enter the unit manually via the input parameters (a2f_unit, dos_unit, Weep_unit, Wen_unit)

## Q2: Why does the order parameter start to oscillate between positive and negative values?
The eliashberg equations are unstable when the coulomb part of the order parameter exceeds the phonon part at any iteration. Double check if the *muc_ME* is set correctly. Try setting the damping of the coulomb part to a higher value via *nItFullCoul* or use a different mixing factor. It is also possible that the material is simply not a superconductor at the given temperature.

## Q3: What does `NaN`, -1 and "" indicate in the input structure?
These sentinel values indicate that either a default value is used or that the input is optional: `NaN` for real-valued fields, `-1` for integer fields (line counts, column indices) and `""` for strings. If such an input is left at its sentinel, it will be determined during the run and the placeholder will be overwritten.

E.g. if the Fermi energy (*ef*) is unspecified (`NaN`), its value will be extracted from the header of the dos-file and the `NaN` will be replaced by the actual value.

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

## Q5: Which role do the cutoffs and grids play in the real-axis solver?
The real-axis solver is far more sensitive to its grids than the imaginary-axis solver: its integrands are complex, sharply peaked and their poles move during the iteration.
Most of the failures reported so far can be traced back to one of the three parameters below.
A detailed description of where each of them enters is given on the [Real Axis Solver](@ref) page.

### The ``\mu``-update diverges — increase *reOmega_c_shift*
The chemical potential is fixed by charge neutrality, and the electron number in the superconducting state is obtained from an ``\omega``-integral of the shift channel.
That integral converges slowly, because ``\chi(\omega)`` decays much more slowly than ``\Delta(\omega)`` or ``Z(\omega)``: this is exactly why ``\chi`` has its own, much larger cutoff (*reOmega_c_shift*, default 25000 meV) instead of sharing *reOmega_c*.

If *reOmega_c_shift* is chosen too small, the tail of that integral is truncated and the computed electron number is systematically wrong.
The root finder then chases a root of a function that never crosses zero within a physical window and walks ``\mu`` away from the Fermi level.
Symptoms are an *ef-mu* column in the console output that grows from iteration to iteration instead of settling, or a hard stop with `muError.png`/`muError.dat` reporting that no root could be found.

**Fix:** increase *reOmega_c_shift* until ``\mu`` settles and the result stops changing. Note that this cutoff is a convergence parameter in its own right and generally has to be much larger than *reOmega_c*.

### ``\mu \to 0`` — decrease *depsilon*
The ``\varepsilon``-integrals are evaluated analytically on a piecewise-linear representation of the DOS, on a uniform grid of step *depsilon*.
The charge-neutrality condition is a difference of two nearly equal electron numbers, so it is only as accurate as the DOS representation around ``\varepsilon_F``.

If the ``\varepsilon``-grid is too sparse, the piecewise-linear DOS smears out the structure near the Fermi level and the difference between the normal and the superconducting electron number collapses.
The root finder then returns ``\mu \approx 0`` regardless of temperature, i.e. the *ef-mu* column stays pinned at (or drifts towards) zero even though the DOS is clearly not constant.

**Fix:** decrease *depsilon* until the shift stabilizes. This is the real-axis counterpart of the *itpStepSize*/*itpBounds* convergence test on the imaginary axis; keep in mind that the cost of an iteration grows with the number of ``\varepsilon``-points.

### No ``T_c`` is found although one is expected — decrease *dOmega*
The kernel is not evaluated pair by pair. Instead, the ``\Omega``-integrals over ``\alpha^2F`` are tabulated on a uniform grid of step *dOmega* and the kernel is reconstructed from these tables by interpolation.
*dOmega* therefore controls how well the phonon structure of ``\alpha^2F`` is resolved in the kernel.

If *dOmega* is too coarse, sharp phonon peaks are smeared out, the effective coupling is underestimated and the gap is suppressed.
The solver then reports ``\Delta(0) <`` *minGap* even at low temperatures and the ``T_c`` search terminates without a transition, although Allen-Dynes (printed at the start of the run) predicts a finite ``T_c``.
An Allen-Dynes ``T_c`` that is orders of magnitude above the ``T_c`` returned by the solver is a good indicator for this.

**Fix:** decrease *dOmega* — in particular for materials with narrow phonon peaks. It is the parameter that dominates the setup cost of each temperature, so it is worth converging it on a single temperature rather than during a full ``T_c`` search.