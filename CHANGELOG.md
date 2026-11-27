# Changelog

## v2.0.0

### Breaking changes

Scripts written for IsoME 1.x keep working unless they touch one of the following. All
three are reported by `arguments()` before the run starts, with the migration step in the
message.

- **`shiftcut` removed.** `encut` is now the single ε-cutoff and bounds `χ` and the electron
  number as well. Drop `shiftcut` and give `encut` the value you used for it. The default of
  `encut` moved from 5000 to 2000 meV accordingly.
- **`-1` is no longer the "not set" sentinel.** `mu`, `muc_AD`, `muc_ME`, `ef`, `efW`,
  `typEl` and `mixing_beta` now use `NaN`. Leave the field out to have it inferred; an
  explicit `-1` is rejected rather than silently used as a μ*, a mixing factor or an energy.
  The integer fields (`nheader_*`, `nfooter_*`, `ind_smear`, `nsmear`, column indices) keep
  `-1`.
- **Field types are concrete.** `temps` is a `Vector{Float64}` (was `Vector{Number}`) and
  `mixing_beta` a `Float64` (was `Number`); the `nheader_*` / `nfooter_*` / `nItFullCoul`
  fields are `Int64`. Integer input converts automatically, so `temps = [5, 10]` is still
  fine. A value that cannot convert is reported by field name.

### Changed

- **`Info.txt` renamed to `Input.txt`.** The file lists the input parameters of the run;
  the new name says so. Scripts that collect `Info.txt` from an output directory need to
  be updated.
- **Output directory.** `outdir` now defaults to `IsoME/` inside the current working
  directory instead of the working directory itself. An existing directory is never written
  into: a run counter is appended, so repeated runs land in `IsoME_1/`, `IsoME_2/`, and so on.
- A unit that cannot be read from a file header is now an error naming the `*_unit` input to
  set. Previously the run stopped to prompt for it on stdin, which hangs in a batch job.
- `mu`, `muc_AD`, `muc_ME`, `mixing_beta`, `conv_thr`, `minGap`, `Tc_tol`, `N_it` and `min_it` are
  rejected when negative, at construction and again at the start of a solve. `NaN` is
  unaffected, so leaving a field out to have it inferred still works.
- Further input checks: `typEl` > 0, `0 < mixing_beta ≤ 1`, valid `temps`, `mu_flag` ∈ {0, 1},
  `spinDos` ∈ {1, 2}. Unsorted `itpBounds` are sorted with a warning. A too small `N_it` is
  raised to `max(min_it, nItFullCoul + 1) + 1`.
- A temperature that reaches `N_it` or becomes `NaN` is reported with a warning, treated as
  normal conducting, and counted in a note in the summary. Convergence also requires a gap
  of at least `minGap`.
- The Fermi energy is read from a header keyword first; a keyword-free header is only
  guessed from, with a warning.
- A positive `ef` (or `efW`) is used as typical electronic energy when `typEl` is unset; the
  value used is printed in the log.

### Added

- **`RealAxisSolver`**, which solves the Eliashberg equations directly on the real frequency
  axis in cDOS+μ, vDOS+μ and vDOS+W. Same `arguments`, same input files as
  `EliashbergSolver`. New inputs: `domega`, `n_cheb`, `depsilon`, `min_it`.
- `flag_acon`: Padé continuation of the converged imaginary-axis solution to real
  frequencies, written to `outdir/ACON/`.
- `a2f_Nef`: rescales α²F to the N(ε_F) of the DOS file when the two were computed with
  different values.
- **`Tc_tol`**: the resolution the Tc search is carried to, in K. The default `1.0` is the
  previous behaviour. 
- A docstring for `arguments`, listing every input group (`?arguments`).

### Fixed

- Allen-Dynes/ML estimates no longer blow up for λ < μ*(1+0.62λ): they are reported as 0 K
  with a warning, and the solver starts from a gap of 10·`minGap`.
- `Wen_unit` is now used for the W energy grid (from the `Wen_file` header, else the unit of W).
- `curve_fit` no longer reports convergence when no step was taken.
- The cDOS seed of the real-axis vDOS solver always uses μ* = 0.12.

## v1.0.5 and earlier

See the release history at <https://github.com/cheil/IsoME.jl/releases>.
