# Real Axis Solver

`RealAxisSolver()` solves the isotropic Eliashberg equations directly on the real frequency axis, following the formulation of [Simon](https://doi.org/10.48550/arXiv.2603.18199).
Unlike the imaginary-axis solver it yields the complex self-energy components ``Z(\omega)``, ``\Delta(\omega)`` and ``\chi(\omega)`` without an analytic continuation, at the price of a considerably more involved numerical treatment: every quantity is complex, the integrands carry poles and branch points close to the real axis, and the kernel has to be rebuilt at every temperature.

The solver supports the same approximations as the imaginary-axis code, selected through `cDOS_flag` and `include_Weep`: cDOS``+\mu``, vDOS``+\mu`` and vDOS``+W``.

## Structure of a calculation
For each temperature the solver performs the following steps:

1. **Kernel setup** (once per temperature). The electron-phonon kernel is assembled from ``\alpha^2F`` by carrying out the ``\Omega``-integration. This is the only step that involves the phonons and by far the most expensive part of the setup.
2. **Self-consistency loop**. Starting from a BCS-like guess (or from the converged solution of the previous temperature), each iteration
   - locates the poles of the ``\omega'``-integrand and, if they moved, rebuilds the pole-adapted part of the ``\omega'``-grid,
   - updates the chemical potential so that charge neutrality is preserved (vDOS only),
   - performs the ``\varepsilon``-integration, which reduces the ``(\omega',\varepsilon)``-integrand to a pure function of ``\omega'``,
   - performs the ``\omega'``-integration against the kernel, giving the new ``Z``, ``\phi`` and ``\chi``,
   - mixes the new and old solutions and checks convergence.
3. **Termination**. A temperature counts as converged once the root-mean-square change of ``\Delta`` relative to ``\Delta(0)`` drops below `conv_thr`. If ``\Delta(0)`` falls below `minGap`, the temperature is counted as normal conducting, which is what brackets ``T_c`` in ``T_c``-search mode.

Throughout the loop the primary variables are ``Z``, ``\chi`` and the order parameter ``\phi``; the gap is derived as ``\Delta = \phi/Z`` after mixing.
This matters in the vDOS``+W`` case, where ``\phi(\omega,\varepsilon)`` is energy dependent and ``\Delta(\omega)`` is defined through ``\phi(\omega,\varepsilon_F)/Z(\omega)``.

The frequency grids on which the self-energy is stored are uniform and fixed for the whole run.
``Z(\omega)`` and ``\Delta(\omega)`` live on `1e-1 : domega : reOmega_c`, whereas ``\chi(\omega)`` uses the same step but the much larger cutoff `reOmega_c_shift`.
Both grids deliberately start at ``0.1`` meV rather than at ``0``: several integrands behave like ``1/\omega``, which makes the result sensitive to the first grid point.

## Kernel
The electron-phonon interaction enters the equations through the kernels ``K^{\pm}(\omega,\omega')`` (supplemental eq. (47) of the reference).
Evaluated naively, each pair ``(\omega,\omega')`` requires a principal-value integral of ``\alpha^2F(\Omega)`` over ``\Omega``, weighted with Bose and Fermi factors — an ``\mathcal{O}(N^2)`` set of integrals for grids of ``N \sim 10^4`` points.
IsoME avoids this in two steps.

### The ``\Omega``-integration
The kernel is not an arbitrary function of two variables.
All ``\Omega``-integrals reduce to two principal-value integrals of the *difference* variable ``x``,
```math
I_1(x) = \mathcal{P}\!\int \! d\Omega \, \frac{G(\Omega)}{\Omega-x}~,\qquad
I_2(x) = \mathcal{P}\!\int \! d\Omega \, \frac{G(\Omega)\, n(\Omega)}{\Omega-x}~,
```
where ``G`` is the interpolated ``\alpha^2F`` and ``n`` the Bose function.
``G`` is nonzero only between the first and the last frequency at which ``\alpha^2F`` exceeds ``10^{-6}``; outside this window the integrands vanish identically, which bounds the ``\Omega``-integration range.

The ``\Omega``-axis is composed of a Chebyshev grid inside a small window around the ``1/\Omega`` singularity, where the integrand varies rapidly and the principal value has to be resolved, and of linear grids outside of it, where the integrand is smooth.
The integration itself is performed with the trapezoidal rule.

``I_1`` and ``I_2`` are tabulated on a uniform ``x``-grid of step `dOmega` that spans ``\pm 2\,\omega_{\max}``, so that every combination ``\pm\omega\pm\omega'`` reachable from the frequency grids falls inside the table.
This is the expensive part of the setup, and it is also the reason the kernel cannot be reused across temperatures: both ``I_2`` and the Fermi factors depend on ``\beta``.

### The linear (difference-grid) representation
From ``I_1`` and ``I_2`` the kernel is assembled into two one-dimensional complex tables ``A(x)`` and ``B(x)``, such that
```math
\mathcal{K}(\omega,\omega') = A(x) + f(\omega')\,B(x)~,\qquad x = \omega'-\omega~,
```
with the physical kernels following by folding ``\omega' \to -\omega'``,
```math
K^{+}(\omega,\omega') = \phantom{-}\mathcal{K}(\omega,\omega') - \mathcal{K}(\omega,-\omega')~,\qquad
K^{-}(\omega,\omega') = -\mathcal{K}(\omega,\omega') - \mathcal{K}(\omega,-\omega')~.
```
Only the ``\mathcal{O}(N)`` tables ``A`` and ``B`` are stored; the ``\mathcal{O}(N^2)`` kernel matrix is never formed.
Evaluating the kernel at an arbitrary pair ``(\omega,\omega')`` then amounts to interpolating ``A`` and ``B`` at ``x = \omega'-\omega`` and ``x = -\omega-\omega'``, plus one Fermi factor.

Because the kernel is reconstructed by *interpolating* these tables, `dOmega` controls how well the ``\Omega``-integral — and with it the phonon structure — is resolved.
A too coarse `dOmega` smears the sharp features of ``\alpha^2F`` and can suppress the gap altogether; see the [FAQ](@ref).

## ``\varepsilon``-integration
In the vDOS approximations the equations contain an integral over the electronic energy ``\varepsilon``, weighted with the density of states ``N(\varepsilon)``.
The ``\varepsilon``-dependence of the integrand is a product of two Lorentzians, centred at ``\varepsilon = -(\chi \pm \varepsilon_p)`` with
```math
\varepsilon_p(\omega') = \sqrt{\omega'^2 Z(\omega')^2 - \phi(\omega')^2}~.
```
Since ``Z``, ``\phi`` and ``\chi`` are complex, these Lorentzians can become arbitrarily narrow wherever their imaginary parts collapse — which is exactly where a numerical quadrature over ``\varepsilon`` would fail.

IsoME therefore evaluates the ``\varepsilon``-integral **analytically**.
The DOS is taken to be piecewise linear between its grid points, ``N(\varepsilon) = M_0 + M_1\,\varepsilon`` on each interval, and the resulting integrals
```math
\int_{\varepsilon_l}^{\varepsilon_r} \frac{\varepsilon^k}{(\varepsilon+a)^2+b^2}\, d\varepsilon~, \qquad k = 0\ldots3
```
have closed-form expressions in terms of `atan` and `log`.
Summing these interval moments, weighted with ``M_0`` and ``M_1``, over all intervals yields the ``\omega'``-integrands ``g_Z``, ``g_\phi`` and ``g_\chi`` exactly.
The only remaining approximation is the piecewise-linear representation of the DOS itself, so the ``\varepsilon``-grid has to resolve the DOS, but not the Lorentzians.

The ``\varepsilon``-grid is uniform with step `depsilon` and spans the range set by `encut` (or the extent of the DOS file, whichever is smaller).
This differs from the imaginary-axis solver, which uses the piecewise interpolation controlled by `itpBounds` and `itpStepSize`.

Two cases are treated separately:

- **cDOS**``+\mu``: there is no ``\varepsilon``-integral at all. The integrands reduce to the closed BCS-like forms ``g_Z = \mathrm{Re}\left[\omega'/\sqrt{\omega'^2-\Delta^2}\right]`` and ``g_\phi = \mathrm{Re}\left[\Delta/\sqrt{\omega'^2-\Delta^2}\right]``.
- **vDOS**``+W``: ``\phi(\omega,\varepsilon)`` is energy dependent and is likewise treated as piecewise linear, ``\phi = \Phi_0 + \Phi_1\varepsilon``. This modifies the Lorentzian centres to
  ```math
  S = \frac{\chi + \Phi_0\Phi_1}{1+\Phi_1^2}~, \qquad
  P = \sqrt{\frac{\omega'^2Z^2 - \chi^2 - \Phi_0^2}{1+\Phi_1^2} + S^2}~,
  ```
  which take over the roles of ``\chi`` and ``\varepsilon_p``. The Coulomb term additionally couples the coefficients of the ``\phi``-branch to the piecewise-linear moments of ``N(\varepsilon')W(\varepsilon,\varepsilon')``, which is evaluated as a matrix product over the ``\varepsilon``-grid.

In all cases the sign of the ``Z``-integrand fixes the branch: it has to be positive for causality, and the ``\phi``- and ``\chi``-integrands are flipped accordingly.

## ``\omega'``-integration
After the ``\varepsilon``-integration the remaining task is
```math
I_Z(\omega) = \int \! d\omega' \, K^{-}(\omega,\omega')\, g_Z(\omega')~, \qquad
I_\phi(\omega) = \int \! d\omega' \, K^{+}(\omega,\omega')\, g_\phi(\omega')~,
```
and analogously for ``\chi`` on the larger frequency grid.
These integrands are sharply peaked: the ``\varepsilon``-integration leaves poles wherever a Lorentzian width collapses, i.e. at the roots of ``\mathrm{Im}(\chi \pm \varepsilon_p)``, plus a branch point at the gap edge, where ``|\varepsilon_p|^2`` becomes minimal.
Their positions move from iteration to iteration as the self-energy changes, so a fixed grid would either be inaccurate or prohibitively dense.

The ``\omega'``-grid is therefore split at ``\omega'_{\max} = 2\Delta(0)`` of the starting state (clamped to a sensible range) into a **head** and a **tail**:

- **Head** — ``(0, \omega'_{\max}]``, the region in which the poles live. Every pole receives a Chebyshev cluster of `n_cheb` points, dense towards the pole from both sides, and the clusters tile the whole interval up to the midpoints between neighbouring poles. In the cDOS case there is a single pole, namely the root of ``\omega' - \mathrm{Re}\,\Delta(\omega')``.
- **Tail** — ``[\omega'_{\max}, \omega_c]``, where the integrand is smooth. A uniform grid of step `domega` up to `reOmega_c` (cDOS) resp. `reOmega_c_shift` (vDOS).

Both parts are integrated with trapezoidal weights, and the junction point belongs to both grids, so the split is exact.

The split is also what makes the solver affordable:

- The head kernel blocks are materialized — they are only as wide as the head is long — and the integration is a matrix-vector product. They are rebuilt only when the poles actually move by more than ``10^{-3}\,\omega'_{\max}`` or when their number changes, not in every iteration. If a pole wanders beyond ``\omega'_{\max}``, the workspace is rebuilt with a larger head.
- The tail shares its step with the ``\omega``-grids. The kernel arguments then depend only on ``j-i`` (difference channel, Toeplitz) and ``i+j`` (sum channel, Hankel), so ``A`` and ``B`` have to be sampled only ``\mathcal{O}(n_\omega + n_{\text{tail}})`` times. The tail contribution is accumulated on the fly from these samples, and its kernel matrix is never built.

## Chemical potential
In the vDOS approximations the transition to the superconducting state shifts the chemical potential.
With `mu_flag = 1` (the default), IsoME determines ``\mu`` in every iteration by requiring that the electron number in the superconducting state matches the one in the normal state.
The root of ``N_e^{\text{nsc}}(\mu) - N_e^{\text{sc}}(\mu)`` is located with a Regula-Falsi method, starting from a ``\pm 50`` meV window around the current Fermi level that is widened in steps of 50 meV until it brackets a sign change.
Only the ``Z``- and the ``\chi``-branch of the ``\varepsilon``-integral are needed for this, so the update is considerably cheaper than a full iteration — in particular, ``W`` never enters charge conservation.

The ``\mu``-update is the part of the real-axis solver that reacts most sensitively to the grids and cutoffs; the [FAQ](@ref) lists the typical failure modes and how to recognize them.

## Self-consistency and convergence
The default mixing is linear with an iteration-dependent factor that ramps from ``1`` down to ``0.5``.
A fixed factor can be enforced through `mixing_beta`, and `broyden_flag = 1` switches to Broyden mixing (with a history depth of `broyden_mem`) of ``Z``, ``\chi`` and ``\phi``.
In the vDOS``+W`` case the Coulomb contribution is ramped up over the first `nItFullCoul` iterations, which stabilizes the early iterations.

vDOS calculations are never started from scratch: the first temperature is seeded with a cDOS solution at the same temperature, and every subsequent temperature starts from the converged solution of the previous one.
``T_c``-search mode uses the same bisection-and-fit strategy as the imaginary-axis solver, so ``T_c`` ends up bracketed by the highest temperature with ``\Delta(0) >`` `minGap` and the lowest one without a solution.

## Relevant input parameters
| Parameter | Role |
|:----------|:-----|
| `reOmega_c`, `domega` | Cutoff and step of the ``Z(\omega)``, ``\Delta(\omega)`` grid, and of the ``\omega'``-tail in cDOS calculations |
| `reOmega_c_shift` | Cutoff of the ``\chi(\omega)`` grid and of the ``\omega'``-tail in vDOS calculations |
| `dOmega` | Step of the kernel tables ``A(x)``, ``B(x)``; sets the resolution of the ``\Omega``-integral |
| `n_cheb` | Chebyshev points per pole in the head region of the ``\omega'``-grid |
| `depsilon` | Step of the ``\varepsilon``-grid (vDOS) |
| `encut` | Range of the ``\varepsilon``-integration |

See [Input](@ref) for defaults and types, and the [FAQ](@ref) for the practical consequences of these choices.

```@docs
RealAxisSolver
```
