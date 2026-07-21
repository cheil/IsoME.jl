# Matsubara Solver

`EliashbergSolver()` solves the isotropic Migdal-Eliashberg equations on the imaginary (Matsubara) frequency axis, as described in the [IsoME](https://doi.org/10.1016/j.cpc.2025.109720) paper.
On the imaginary axis all self-energy components are real and smooth, which makes this solver fast and robust — the natural choice for ``T_c`` searches and high-throughput work.
The price is that spectral quantities are not obtained directly; they require an analytic continuation to real frequencies.

The approximation is selected through `cDOS_flag` and `include_Weep`: cDOS``+\mu``, vDOS``+\mu``, vDOS``+W`` (and, for completeness, cDOS``+W``, which is not recommended).

## Structure of a calculation
For each temperature the solver

1. builds the Matsubara grid and precomputes the electron-phonon coupling ``\lambda(i\omega_n - i\omega_m)``,
2. initializes ``\Delta``, ``Z`` (and ``\chi``, ``\phi_c`` where applicable) from a BCS-like guess,
3. iterates the Eliashberg equations to self-consistency, updating the chemical potential and mixing the solutions on the way,
4. optionally continues the converged solution to the real axis.

The self-energy components are ``Z(i\omega_n)``, ``\Delta(i\omega_n)`` and, in vDOS mode, the shift ``\chi(i\omega_n)``.
In the ``W`` approximations the order parameter is split into a phononic part ``\phi_{ph}(i\omega_n)`` and a Coulomb part ``\phi_c(\varepsilon)``, and the gap follows as ``\Delta = (\phi_{ph}+\phi_c)/Z``.
The equations follow Pickett, PRB **26**, 1186 (1982) and Lee *et al.*, npj Comput. Mater. **9**, 156 (2023) for the full-bandwidth case, and Margine & Giustino, PRB **87**, 024505 (2013) for the constant-DOS case.

## Matsubara grid
The fermionic Matsubara frequencies ``\omega_n = (2n+1)\pi k_B T`` are generated up to the cutoff `imOmega_c`, which fixes their number as
```math
M = \left\lceil \frac{1}{2}\left(\frac{\omega_c}{\pi k_B T} - 1\right) \right\rceil~.
```
The grid is thus temperature dependent: the lower the temperature, the denser the frequencies and the more of them fall below the cutoff.
This is what makes low-temperature calculations expensive, and it is the reason for sparse sampling (see below).

`imOmega_c` is a convergence parameter, but it is bounded from above in practice: ``\mu^*_{ME}`` is adapted to the cutoff through the Morel-Anderson formula, and that adaptation breaks down for very large cutoffs.
See [Best practices](@ref).

The electron-phonon coupling
```math
\lambda(i\omega_n - i\omega_m) = \int \! d\Omega \, \frac{2\,\Omega\,\alpha^2F(\Omega)}{\Omega^2 + (\omega_n-\omega_m)^2}
```
depends only on the difference of two Matsubara indices. It is therefore precomputed once per temperature for all ``2M+2`` distinct values by trapezoidal integration over the ``\alpha^2F`` grid, and the Matsubara sums then reduce to dot products with index-shifted slices of this table.

## Energy grid
In the vDOS and ``W`` approximations the equations contain an integral over the electronic energy ``\varepsilon``, weighted with the DOS.
On the imaginary axis the integrand is a smooth Lorentzian-like function of ``\varepsilon``, so the integral is done numerically with the trapezoidal rule — no analytic treatment is needed, in contrast to the real-axis solver.

The accuracy is set entirely by the ``\varepsilon``-grid, which is built by interpolating the DOS (and ``W``) onto a piecewise-uniform mesh:
`itpBounds` defines regions around the Fermi level, and `itpStepSize` gives the step used within each of them.
With the defaults (`itpBounds = [100, 500]`, `itpStepSize = [1, 5, 50]`) the grid uses a 1 meV step within ``\pm 100`` meV of ``\varepsilon_F``, 5 meV out to ``\pm 500`` meV, and 50 meV beyond that.
The structure of the DOS near ``\varepsilon_F`` is what drives the results, hence the fine step there and the coarse one in the wings.

Two cutoffs limit the ``\varepsilon``-integrations:

- `encut` bounds the grid itself, and with it the ``Z``- and ``\Delta``-integrations. It is a scalar giving the symmetric window ``\pm`` `encut`.
- `shiftcut` bounds the ``\chi``-integration and the charge-neutrality condition, and should be chosen smaller than `encut`. The shift converges faster in ``\varepsilon`` than the other components, and the electron number is only meaningful within a window in which the DOS is well resolved.

## Sparse sampling
As ``T \to 0`` the number of Matsubara frequencies below the cutoff grows like ``1/T``, while the self-energy components stay smooth functions of ``\omega_n``.
Evaluating the equations at every frequency is therefore wasteful at low temperatures.

For ``T <`` `sparseSamplingTemp` IsoME evaluates the equations only on a sparse subset of frequencies and reconstructs the rest by linear interpolation.
The subset is the Matsubara sampling of the intermediate representation (IR) basis, obtained from a `FiniteTempBasis` for the given ``\beta`` and cutoff, and restricted to the frequencies below ``M+1`` (with the last frequency added explicitly, so the interpolation never has to extrapolate towards the cutoff).
Only the outer loop over ``\omega_n`` is sparsified: the Matsubara sums and ``\varepsilon``-integrals inside each equation still run over the full grid, so the sparse sampling reduces the cost without changing the physics.
The number of sampled frequencies is reported in the log file next to the total.

## Chemical potential
In the vDOS approximations the chemical potential is updated in every iteration (from the second one on) when `mu_flag = 1`, by requiring the electron number to be the same in the superconducting and the normal state.
The electron number in the superconducting state is a Matsubara sum over the DOS-weighted Green's function; the resulting root of ``N_e^{\text{nsc}}(\mu) - N_e^{\text{sc}}(\mu)`` is found with a Regula-Falsi method, starting from a ``\pm 50`` meV window around the current Fermi level that is widened until it brackets a sign change.
Both electron numbers are integrated over the `shiftcut` window rather than the full grid.

If the electron number *increases* with decreasing ``\mu``, or if no sign change is found within ``\pm 5`` eV, the run stops and writes `muError.png`/`muError.dat` for inspection.
This almost always points at the DOS input file rather than at the solver.

## Self-consistency and convergence
The mixing is linear with an iteration-dependent factor ramping from ``1`` down to ``0.5``, unless `mixing_beta` fixes it.
The Coulomb contribution is ramped up over the first `nItFullCoul` iterations, which prevents the instability that arises when the Coulomb part of the order parameter exceeds the phononic part early in the iteration (see the [FAQ](@ref)).

Convergence is measured by the relative change of the gap, ``\sum_n |\Delta_n^{i} - \Delta_n^{i-1}| / \sum_n |\Delta_n^{i}|`` (at ``\varepsilon_F`` in the ``W`` approximations), and is accepted once it drops below `conv_thr` and at least `min_it` iterations have passed.
A temperature is considered normal conducting when ``\Delta(i\omega_0)`` falls below `minGap`, which is what brackets ``T_c`` in ``T_c``-search mode: the search starts from a machine-learning estimate, brackets the transition by bisection, and accelerates convergence with a fit of ``\Delta(T)``.

## Analytic continuation
With `flag_acon = true` the converged imaginary-axis solution is continued to real frequencies with Padé approximants.
``\Delta``, ``Z`` and ``\chi`` are continued separately (which is more stable than continuing the Green's function directly), and the resulting real-frequency gap, renormalization and quasiparticle DOS are plotted and written out.
Padé is sensitive to the number of input frequencies and to noise in the converged solution; if the continued quantities look unphysical, the direct [Real Axis Solver](@ref) is the more reliable route.

## Relevant input parameters
| Parameter | Role |
|:----------|:-----|
| `imOmega_c` | Matsubara cutoff; sets the number of frequencies at a given temperature |
| `encut` | Range of the ``\varepsilon``-grid and of the ``Z``-/``\Delta``-integrations |
| `shiftcut` | Range of the ``\chi``- and charge-neutrality integrations |
| `itpBounds`, `itpStepSize` | Regions and step sizes of the DOS/``W`` interpolation around ``\varepsilon_F`` |
| `sparseSamplingTemp` | Temperature below which IR sparse sampling is used |
| `flag_acon` | Analytic continuation of the converged solution |

See [Input](@ref) for defaults and types, and [Best practices](@ref) for guidance on convergence tests.

```@docs
EliashbergSolver
```
