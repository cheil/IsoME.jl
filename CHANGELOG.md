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

- **Output directory.** `outdir` now defaults to `IsoME/` inside the current working
  directory instead of the working directory itself. An existing directory is never written
  into: a run counter is appended, so repeated runs land in `IsoME_1/`, `IsoME_2/`, and so on.
- A unit that cannot be read from a file header is now an error naming the `*_unit` input to
  set. Previously the run stopped to prompt for it on stdin, which hangs in a batch job.
- `mu`, `muc_AD`, `muc_ME`, `mixing_beta`, `conv_thr`, `minGap`, `N_it` and `min_it` are
  rejected when negative, at construction and again at the start of a solve. `NaN` is
  unaffected, so leaving a field out to have it inferred still works.

### Added

- **`RealAxisSolver`**, which solves the Eliashberg equations directly on the real frequency
  axis in cDOS+μ, vDOS+μ and vDOS+W. Same `arguments`, same input files as
  `EliashbergSolver`. New inputs: `domega`, `n_cheb`, `depsilon`, `min_it`.
- `min_it` is honoured by both solvers. The imaginary-axis solver previously used a fixed
  minimum of 15 iterations, so with the default `min_it = 10` it may now converge earlier.
- `flag_acon`: Padé continuation of the converged imaginary-axis solution to real
  frequencies, written to `outdir/ACON/`.
- `a2f_Nef`: rescales α²F to the N(ε_F) of the DOS file when the two were computed with
  different values.
- A docstring for `arguments`, listing every input group (`?arguments`).

## v1.0.5 and earlier

See the release history at <https://github.com/cheil/IsoME.jl/releases>.
