# Troubleshooting

This page explains the error and warning messages IsoME can print, what causes them, and how to
react. If you only have a conceptual question, the [FAQ](@ref) is the better starting point;
this page is meant to be searched for the exact text of a message you are seeing.

## How IsoME reports problems

IsoME distinguishes two kinds of messages:

- **Fatal errors** stop the run. They are reported as `<stage>. Stopping now!`, where the stage
  tells you *where* it happened — `in input structure`, `while reading the inputs`, or
  `while solving the Eliashberg equations`. The full backtrace is written to a `CRASH` file in the
  output directory; the console and log only show the message. Start from the message below, and
  open the `CRASH` file if you need the exact call stack.
- **Warnings** are printed but the run continues. They flag a fallback (a default value was
  substituted) or an optional step that failed (a plot, an output file). The final result is still
  produced.

Some checks also drop diagnostic artifacts into the output directory: `muError.png` / `muError.dat`
and `muError_shift.png` / `muError_shift.dat` for a failed [μ-update](@ref "The μ-update"), and
`muc_ME.png` for the [μ\* conversion](@ref "μ* conversion").

### The output directory contains only a log.txt

No `CRASH` file and no `Summary.dat` — this looks like a silent failure, but it is not: it is what
the output directory looks like whenever a run *ends without reaching its last step*. The three
files are written at different times, which is what makes the state readable:

| File | Written |
|------|---------|
| `log.txt` | opened when the run starts, appended to throughout |
| `CRASH` | whenever an exception was caught — a fatal error, or an optional step that failed (see [Optional output steps failed](@ref)) |
| `Input.txt`, `Summary.dat`, figures | only after the equations are solved, at the very end of a successful run |

So a lone `log.txt` means the run neither finished nor hit a handled error. The log itself cannot
tell you which: a run that is still working and a run that was stopped leave the same partially
written file, both ending in the middle of a temperature step. **Check whether the calculation is
still running** — `ps`/`htop` for a local run, `condor_q` or the equivalent for a batch job:

- **It is still running.** Nothing is wrong. `Input.txt`, `Summary.dat` and the figures are written
  in one go after the equations are solved, so their absence says nothing about progress while the
  temperature loop is still printing. This is the normal state of every output folder of a batch
  that has not drained yet.
- **It is no longer running.** The run was stopped before it reached its last step. There are two
  common reasons, neither of which produces a `CRASH` file:
  - **Interrupted by the user (`Ctrl+C`).** An interrupt is not treated by IsoME's error handling —
    it is an abort, not a failure, so no crash is reported and nothing is written. Simply start the
    run again; there is nothing to fix in the input.
  - **Killed from the outside.** The OOM killer, a batch system evicting or holding the job, an
    explicit `kill`, or the machine going down. Julia never gets to run any handler in that case
    either. Check the scheduler's own files (for HTCondor the `.log` / `.err` of the job) for the
    reason, and resubmit — with more memory if the job was killed for its size.

The same reasoning applies to a partially written output directory (`Input.txt` present,
`Summary.dat` missing): everything up to the point of the interruption was written, everything
after it was not.

## Input and file-reading errors

These are raised `while reading the inputs`. They almost always mean a path, a flag, or a file
header needs to be corrected — see the [Input](@ref) page for the parameters and expected file
formats, and the [FAQ](@ref) for the assumed column layout.

| Message | Cause | What to do |
|---------|-------|------------|
| `Invalid path to a2f-file!` | `a2f_file` does not point to an existing file | Check the path; use an absolute path if in doubt. |
| `Invalid path to Dos-file!` | `dos_file` missing, but required (`cDOS_flag = 0` or `include_Weep = 1`) | Provide a DOS file, or switch to `cDOS_flag = 1` if a constant DOS is intended. |
| `Invalid path to Weep or Wen-file!` | `include_Weep = 1` but `Weep_file` (or `Wen_file`) is missing | Provide the file, or set `include_Weep = 0` to use the μ approximation. |
| `Invalid cDOS_flag value. …` | `cDOS_flag` is not 0 or 1 | Use `0` (vDOS) or `1` (cDOS). |
| `Invalid include_Weep value. …` | `include_Weep` is not 0 or 1 | Use `0` (μ approximation) or `1` (W(ε, ε′)). |
| `The real-axis solver only supports the vDOS+W approximation …` | `RealAxisSolver` called with `include_Weep = 1` and `cDOS_flag = 1` | Set `cDOS_flag = 0`. The real-axis solver has no cDOS+W mode. |
| `Unknown mode! …` | `cDOS_flag` / `include_Weep` combination is not a supported mode | Check both flags against the modes in the [Input](@ref) page. |
| `Couldn't write into <outdir>! …` | `outdir` is not writable, or the path is invalid | Fix the path or the permissions. Missing parent directories are created automatically. |
| `Could not determine the unit of the <file>-file. …` | The unit in the file header was not recognized | Units are case-sensitive; supported are **meV, eV, THz, Ry, Ha**. Set it manually via the matching `*_unit` flag (`a2f_unit`, `dos_unit`, `Weep_unit`, `Wen_unit`). See [FAQ Q1](@ref "FAQ"). |
| `Error while reading the fermi energy from the <file>-file.` | The header line for the Fermi energy could not be parsed | Set it manually via `ef` / `efW`, or check the file header. |
| `Could not extract the Fermi energy from the <file>-file.` | Parsing ran but matched no Fermi-energy value | Same as above. |
| `The first column of the <file>-file holds no numeric entry …` | Header/footer auto-detection found no numbers in the first column | Check the column layout of the file, or set the header/footer size manually via `nheader_*` / `nfooter_*`. |
| `α²F stays below 1e-2 over the whole frequency range …` | The selected ``\alpha^2F`` column is (almost) zero everywhere | Check the ``\alpha^2F`` file and the selected smearing column `ind_smear`; the unit conversion may also have scaled the values away. |
| `The a2F frequency grid does not reach 2·domega = …` | `domega` is larger than half the highest frequency of the ``\alpha^2F`` grid | Reduce `domega`, or check that the ``\alpha^2F`` file and its unit are read correctly. |
| `The density of states is zero over the whole energy range.` | The column read as the DOS contains only zeros | Check `dos_file`, its column layout and `spinDos`. |
| `a2f_Nef must be positive (got …).` | `a2f_Nef` was set to a non-positive value | Pass the ``N(\varepsilon_F)`` used for the ``\alpha^2F`` calculation, or leave it at `NaN` to disable the rescaling. |

### Errors from `arguments()` itself

The keywords of `arguments` are checked against the field names and the fields are concretely typed,
so a wrong input is caught when the input structure is built — before any solver runs. Four messages
come from there:

- *unknown input to `arguments`* — the keyword is not a field. Every offending name is listed with a
  "did you mean …?" suggestion, which is what a renamed or misspelled field looks like. An input that
  was removed in 2.0 (e.g. `shiftcut`) instead gets a line naming what replaces it.
- *invalid input to `arguments`* — a required keyword (`a2f_file`) is missing, or a value cannot be
  converted to the declared type. All offending fields are reported at once, each with its expected
  type and a hint (e.g. *pass a vector, e.g. temps = [10.0]*).
- *invalid input to `arguments`*, reporting a negative value — `mu`, `muc_AD`, `muc_ME`,
  `mixing_beta`, `conv_thr`, `minGap`, `N_it` and `min_it` have no meaning below zero and are
  rejected here, and again at the start of a solve in case the field was assigned afterwards.
- *outdated input value* — an explicit `-1` on `mu`, `muc_AD`, `muc_ME`, `ef`, `efW`, `typEl` or
  `mixing_beta`. Before IsoME 2.0 these used `-1` as the "not set" sentinel; it is now `NaN`, so a
  `-1` would be taken at face value. Leave the field out, or pass `NaN`.

`?arguments` lists every field with its type.

## The μ-update

In the vDOS+μ and vDOS+W approximations the chemical potential ``\mu`` is fixed at every temperature
by charge neutrality: the electron number in the superconducting state must match the normal state.
IsoME solves this by finding the root of ``N_e^{\text{nsc}}(\mu) - N_e^{\text{sc}}(\mu)``. Two fatal
errors come from this root find, both raised `while solving the Eliashberg equations`:

| Message | Meaning |
|---------|---------|
| `The number of electrons decreases with increasing mu!` | ``N_e(\mu)`` is not monotonically increasing, so the root find has no well-defined bracket. |
| `Error in mu update - Couldn't find a root in the interval [.,.].` | No sign change was found within the searched ``\mu`` window. |

Both write two diagnostics to the output directory:

- `muError.png` / `muError.dat` — ``N_e^{\text{nsc}} - N_e^{\text{sc}}`` against ``\mu``. Inspecting
  this curve is the fastest way to see whether a root exists at all and whether the slope has the
  expected sign.
- `muError_shift.png` / `muError_shift.dat` — the shift channel ``\chi(\omega)`` over its frequency
  grid. ``N_e`` is obtained from an ``\omega``-integral of ``\chi``, so if ``\chi`` has **not decayed
  to ``\approx 0`` at the edge of its grid** the integral is truncated and the update cannot
  converge. This is the direct symptom of a too-small cutoff (real axis: increase
  `omega_c`, see below).

On the **real axis**, the electron number is obtained from an ``\omega``-integral of the shift
channel, and the accuracy of the μ-update is governed by two convergence parameters:

- **The μ-update diverges — increase `omega_c`.**
  The shift ``\chi(\omega)`` decays more slowly than ``\Delta(\omega)`` or ``Z(\omega)``, so it is
  usually ``\chi`` that dictates how large the cutoff has to be. If `omega_c` is too small, the
  tail of the integral is truncated, the electron number is systematically wrong, and the root
  finder walks ``\mu`` away from the Fermi level. The *ef-mu* column in the console then grows from
  iteration to iteration instead of settling.
  **Fix:** increase `omega_c` until ``\mu`` settles, and check `muError_shift.png` to confirm that
  ``\chi`` has decayed to ``\approx 0`` at the edge of the grid. ``Z``, ``\Delta`` and ``\chi`` share
  the single cutoff `omega_c`, so there is no separate knob for the shift channel.

- **``\mu \to 0`` — decrease `depsilon`.**
  The ``\varepsilon``-integrals are evaluated analytically on a piecewise-linear DOS on a uniform
  grid of step `depsilon`. Charge neutrality is a difference of two nearly equal electron numbers,
  so it is only as accurate as the DOS near ``\varepsilon_F``. If the grid is too sparse, that
  difference collapses and the root finder returns ``\mu \approx 0`` regardless of temperature.
  **Fix:** decrease `depsilon` until the shift stabilizes (this is the real-axis counterpart of the
  `itpStepSize` / `itpBounds` convergence test on the imaginary axis).

On the **imaginary axis** the same charge-neutrality root find is solved, but from a Matsubara sum
rather than an ``\omega``-integral. Its convergence likewise depends on the Matsubara cutoff and on
the ``\varepsilon``-grid (the interpolated DOS grid set by `itpBounds` / `itpStepSize`): too small a
cutoff or too coarse a grid can prevent the μ-update from settling. The `muError_shift.png` plot
(here ``\chi`` against the Matsubara frequencies) is written in the same way and is the first thing
to check.

<!-- Imag Axis: f(mu) oscillates "around" monotonic function epsilon grid too coarse; omega_c vs encut unclear -->

### ``\omega`` must be sufficiently larger than ``\varepsilon`` (vDOS)

The μ-update and the shift ``\chi(\omega)`` are built from ``\omega'``-integrals over the
``\varepsilon``-resolved kernels. Their correct limiting behaviour — the ``\omega'``-integral
``J_\omega`` tending to ``1`` for large ``\varepsilon`` — is only recovered in the limit
``\omega \to \infty``, i.e. as long as ``\omega`` is *sufficiently larger* than ``\varepsilon``. If
the ``\varepsilon``-window instead reaches up to (or beyond) the ``\omega``-cutoff, ``J_\omega``
crosses over towards ``0`` for the outermost energies, and both the electron number entering the
μ-update and ``\chi`` acquire a systematic error there.

The crossover is gradual, so no sharp ratio ``\omega/\varepsilon`` can be derived at which
``J_\omega`` goes from ``1`` to ``0``. IsoME therefore requires that the ``\varepsilon``-window stay
inside the frequency cutoff: in vDOS calculations `encut` must not exceed `omega_c`, on the
imaginary axis as well as on the real one. If it does, the run prints

```
encut = … exceeds omega_c = …; the ε-integration has to stay inside the frequency cutoff.
Setting encut = …
```

and continues with the reduced ``\varepsilon``-window. If you need the larger `encut`, do **not**
work around the warning by shrinking the energy window — raise `omega_c` instead (a larger cutoff
is beneficial for the μ-update anyway, see above). Note that ``J_\omega`` is only fully recovered
for ``\omega \gg \varepsilon``, so `encut` close to `omega_c` still leaves a systematic error at
the outermost energies; treat `encut` as a convergence parameter rather than pushing it to the
limit.


<!-- TODO(user): the exact imaginary-axis failure modes (Matsubara cutoff vs ε-grid) still need to
     be pinned down from further testing before this note can be made more prescriptive. -->
<!-- TODO(user): additional μ-update guidance to be provided — see open question. -->
<!-- encut only sensible up to moderate values (2 eV?: Too much energy dependence is neglected in the equations. E.g. a2F is assumed to be constant wrt epsilon. If chi gets huge, most certainly the encut is too large.  -->
<!-- Give LaBH8 as example where larger encut makes results worse -> because of a2F =/= a2F(epsilon) -->

<!-- For all convergence parameter: domega, omega_c, encut, depsilon give examples when they are diverging.  -->

## μ\* conversion

IsoME works with three Coulomb parameters: ``\mu``, and the Morel-Anderson pseudopotentials
``\mu^*_{AD}`` (for Allen-Dynes estimates) and ``\mu^*_{ME}`` (for the Migdal-Eliashberg solver).
When some of them are left at their `NaN` default, the missing ones are derived from the others via
the relations documented on the [Input](@ref) page (Pseudopotentials section). The following
warnings flag a fallback in that conversion — the run continues with the substituted value:

| Warning | Meaning | React by |
|---------|---------|----------|
| `Unable to calculate μ* from μ without a typical electron energy!` | ``\mu \to \mu^*_{AD}`` needs a typical electronic energy, but none was available | Set `typEl` (or `ef` / `efW`). Otherwise `μ*_AD = 0.12` is used. |
| `Couldn't calculate a reasonable μ*_ME from μ*_AD.` | The derived ``\mu^*_{ME}`` fell outside a physical range | Check `muc_ME.png`; set ``\mu^*`` manually or change the Matsubara cutoff. `μ*_ME = min(3·μ*_AD, 0.8)` is used. |
| `Matsubara cutoff would lead to μ*_ME > 4*μ.` | The requested cutoff drives ``\mu^*_{ME}`` above ``4\mu`` | A smaller cutoff is used for the μ → ``\mu^*_{ME}`` conversion only; `omega_c` itself is unchanged and the solver still runs at it. Check `muc_ME.png` and the typical electronic energy `typEl`. |

`muc_ME.png` shows ``\mu^*_{ME}`` as a function of the Matsubara cutoff, which is the quickest way to
judge whether the substituted value is sensible.

One **fatal** error comes from the same conversion:

| Error | Meaning | React by |
|-------|---------|----------|
| `The conversion from μ to μ* gave a negative pseudopotential:` | The Morel-Anderson denominator ``1+\mu\ln(\varepsilon_{el}/\omega)`` has passed through zero, so the resulting ``\mu^*`` is meaningless and would act as an *attractive* Coulomb interaction | Check `typEl` first — it has to be in meV, and a value left in eV is the usual cause — then `omega_c` and `mu`. Setting `muc_AD` / `muc_ME` directly skips the conversion. |

<!-- TODO(user): additional μ* conversion guidance to be provided — see open question. -->

## Real-axis grids and the kernel

The real-axis solver is far more sensitive to its grids than the imaginary-axis solver: its
integrands are complex, sharply peaked, and their poles move during the iteration. Besides the
μ-update parameters above, the ``\omega``-grid as a whole — its cutoff `omega_c` and its step
`domega` — can be at fault:

- **No ``\mathrm{T}_C`` is found although one is expected — check the ``\omega``-grid.**
  The gap comes out below `minGap` even at low temperature and the ``\mathrm{T}_C`` search ends without a
  transition, although the Allen-Dynes ``\mathrm{T}_C`` printed at the start of the run is finite; a large
  gap between the two is the indicator. Too small an `omega_c` truncates the ``\omega'``-integration,
  and too coarse a `domega` smears the sharp features of ``\alpha^2F`` and underestimates the
  coupling, since the kernel tables share the `domega` lattice with the ``\omega``- and
  ``\omega'``-grids. **Fix:** raise `omega_c` and/or decrease  `domega`. Building the kernel dominates the setup cost of each temperature, so converge on a
  single temperature rather than during a full ``\mathrm{T}_C`` search.

A fuller description of where each of these parameters enters is on the
[Real Axis Solver](@ref) page.

## The interpolation grid (imaginary axis)

The imaginary-axis solver interpolates the DOS (and ``W``) onto a piecewise-uniform energy grid
before solving, controlled by `itpBounds` and `itpStepSize`. `itpBounds` lists the boundaries of the
regions around the Fermi level, and `itpStepSize` the step used in each region — so `itpStepSize`
must have **exactly one more entry than** `itpBounds`: one step per region, plus one for the range
beyond the outermost bound.

| Message | Cause | Fix |
|---------|-------|-----|
| `Number of interpolation steps (itpStepSize) and interpolation bounds (itpBounds) do not match` | `length(itpStepSize) != length(itpBounds) + 1` | Add or remove a step so `itpStepSize` has one more entry than `itpBounds`. |

With the defaults, the grid uses a 1 meV step within ``\pm 100`` meV of ``\varepsilon_F``, 5 meV up
to ``\pm 500`` meV, and 50 meV out to `encut` (three step sizes, two bounds). See the
[Input](@ref) page (Interpolation of the energy grid) for the convergence test.

## Optional output steps failed

Messages such as *Error while printing the summary.*, *Error while creating the Info file.*,
*Error while plotting. Skipping plots.*, *Error while saving self energy components.*, or
*Error in analytic continuation.* are **warnings**, not failures. These steps run after the physics
is already solved, so the numerical result is unaffected and the run continues. If one of them keeps
failing across runs, it points to a bug or an environment problem (e.g. a plotting backend): please
open an issue at <https://github.com/cheil/IsoME.jl/issues> and attach the `CRASH` file.

## Internal errors

A handful of messages guard internal invariants — grid junctions, kernel block shapes, root
brackets in the mixing routine (`No real root in [a,b]`, `max number of iterations exceeded`,
`head and tail grids must share the junction point`, dimension mismatches, and similar). These
should never be triggered by valid input; they indicate a bug rather than something to fix in your
setup. If you hit one, please open an issue at <https://github.com/cheil/IsoME.jl/issues> and attach
the `CRASH` file together with the input that produced it.
