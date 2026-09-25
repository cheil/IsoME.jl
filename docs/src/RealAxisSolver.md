# Real Axis Solver

`RealAxisSolver()` solves the isotropic Eliashberg equations directly on the real frequency axis, following the formulation of [Simon *et al.*](https://doi.org/10.48550/arXiv.2603.18199).
Unlike the imaginary-axis solver it yields the complex self-energy components ``Z(\omega)``, ``\Delta(\omega)`` and ``\chi(\omega)`` without an analytic continuation.
That comes at the cost of a considerably more involved numerical treatment: every quantity is complex, the integrands carry poles and branch points close to the real axis, and the kernel has to be rebuilt at every temperature.

The solver supports the same approximations as the imaginary-axis code, selected through `cDOS_flag` and `include_Weep`: cDOS``+\mu``, vDOS``+\mu`` and vDOS``+W``.

## Structure of a calculation
For each temperature the solver goes through three steps:

1. **Kernel setup** (once per temperature). The electron-phonon kernel is assembled from ``\alpha^2F`` by carrying out the ``\Omega``-integration. This is the only step that involves the phonons and by far the most expensive part of the setup.
2. **Self-consistency loop**. Starting from a BCS-like guess (or from the converged solution of the previous temperature), each iteration
   - locates the poles of the ``\omega'``-integrand and, if they moved, rebuilds the pole-adapted part of the ``\omega'``-grid,
   - updates the chemical potential so that charge neutrality is preserved (vDOS only),
   - performs the ``\varepsilon``-integration, which reduces the ``(\omega',\varepsilon)``-integrand to a pure function of ``\omega'``,
   - performs the ``\omega'``-integration against the kernel, giving the new ``Z``, ``\phi`` and ``\chi``,
   - mixes the new and old solutions and checks convergence.
3. **Termination**. A temperature counts as converged once the relative change of the gap, ``\sum_\omega |\Delta^{i}-\Delta^{i-1}| / \sum_\omega |\Delta^{i}|``, drops below `conv_thr` and at least `min_it` iterations — and one more than `nItFullCoul` — have passed. If the gap falls below `minGap`, the temperature is counted as normal conducting, which is what brackets ``\mathrm{T}_C`` in ``\mathrm{T}_C``-search mode.

Throughout the loop the primary variables are ``Z``, ``\chi`` and the order parameter ``\phi``; the gap is derived as ``\Delta = \phi/Z`` after mixing.
This matters in the vDOS``+W`` case, where ``\phi`` is split into a phononic part ``\phi_{ph}(\omega)`` and an energy-dependent Coulomb part ``\phi_c(\varepsilon)``, and the gap follows as ``\Delta(\omega) = [\phi_{ph}(\omega)+\phi_c(\varepsilon_F)]/Z(\omega)``.

The frequency grids on which the self-energy is stored are uniform and fixed for the whole run.
``Z(\omega)``, ``\Delta(\omega)`` and ``\chi(\omega)`` all live on the same grid of step `domega` up to `omega_c`, and the ``\omega'``-integration covers the same range.
The grid deliberately does not start at ``0``, but at the smallest multiple of `domega` that is at least ``0.1`` meV: several integrands behave like ``1/\omega``, which makes the result sensitive to the first grid point.
The gap reported in the console, in `Summary.dat` and to the ``\mathrm{T}_C`` search is read off at the **gap edge** ``\omega_g`` — the root of ``\omega - \mathrm{Re}\,\Delta(\omega)``.

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
The support of ``G`` reaches up to the last frequency at which ``\alpha^2F`` exceeds ``10^{-2}``, and starts at the smaller of the first such frequency and the first grid point above ``2\,``​`domega`; outside this window the integrands vanish identically, which bounds the ``\Omega``-integration range.

The ``\Omega``-axis is composed of a Chebyshev grid inside a ``\pm 3`` meV window around the ``1/\Omega`` singularity, where the integrand varies rapidly and the principal value has to be resolved, and of linear grids out to the width over which ``G`` is nonzero, where the integrand is smooth (300 points each).
The integration itself is performed with the trapezoidal rule.
Only the few ``x`` whose pole falls inside the support of ``G`` need a point-by-point evaluation; for all others the ``\Omega``-dependent part factorizes and each ``x`` reduces to a dot product.

``I_1`` and ``I_2`` are tabulated on a uniform ``x``-grid of step `domega` that spans ``\pm 2\,\omega_{\max}``, so that every combination ``\pm\omega\pm\omega'`` reachable from the frequency grids falls inside the table.
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

The tables share their step with the ``\omega``- and ``\omega'``-grids: they are only ever *sampled* at multiples of `domega` — the tail is a `domega` grid and the head lattice is ``domega\cdot\mathbb{Z}`` — so a finer kernel step would only refine a table that is subsampled again afterwards.
`domega` therefore also controls how well the ``\Omega``-integral, and with it the phonon structure, is resolved.
A too coarse `domega` smears the sharp features of ``\alpha^2F`` and can suppress the gap altogether; see the [FAQ](@ref).

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
In vDOS runs `encut` has to stay inside `omega_c`: a larger value is clamped down to it with a warning, so raise `omega_c` if you need a wider ``\varepsilon``-window.
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
and analogously for ``\chi``, on the same frequency grid.
These integrands are sharply peaked: the ``\varepsilon``-integration leaves poles wherever a Lorentzian width collapses, i.e. at the roots of ``\mathrm{Im}(\chi \pm \varepsilon_p)``, plus a branch point at the gap edge, where ``|\varepsilon_p|^2`` becomes minimal.
Their positions move from iteration to iteration as the self-energy changes, so a fixed grid would either be inaccurate or prohibitively dense.

The ``\omega'``-grid is therefore split at ``\omega'_{\max} = 2\,\omega_g`` of the starting state (clamped to ``[10\,`` `domega` ``,\ `` `omega_c` ``/2]``) into a **head** and a **tail**:

- **Head** — ``(0, \omega'_{\max}]``, the region in which the poles live. Every pole receives a Chebyshev cluster of `n_cheb` points, dense towards the pole from both sides, and the clusters tile the whole interval up to the midpoints between neighbouring poles. In the cDOS case there is a single pole, namely the root of ``\omega' - \mathrm{Re}\,\Delta(\omega')``.
- **Tail** — ``[\omega'_{\max}, \omega_c]``, where the integrand is smooth. A uniform grid of step `domega` up to `omega_c`, in both cDOS and vDOS.

Both parts are integrated with trapezoidal weights, and the junction point belongs to both grids, so the split is exact.

The split is also what makes the solver affordable:

- The head is contracted through a *hat lattice*: because ``A`` and ``B`` are piecewise linear on the ``domega``-lattice, the sum over the Chebyshev nodes can be scattered onto ``M \approx \omega'_{\max}/`` `domega` lattice cells first, after which the kernel is read straight off the tables at exact multiples of `domega` — no interpolation and no kernel block. The lattice depends only on where the head nodes sit, so it is rebuilt (in ``\mathcal{O}(n_{\text{head}})``, without touching any kernel data) only when the poles actually move by more than ``10^{-3}\,\omega'_{\max}`` or when their number changes. If a pole wanders beyond ``\omega'_{\max}``, the workspace is rebuilt with a larger head.
- The tail shares its step with the ``\omega``-grids. The kernel arguments then depend only on ``j-i`` (difference channel, Toeplitz) and ``i+j`` (sum channel, Hankel), so ``A`` and ``B`` have to be sampled only ``\mathcal{O}(n_\omega + n_{\text{tail}})`` times. The tail contribution is accumulated on the fly from these samples, and its kernel matrix is never built.

Since ``Z``, ``\phi`` and ``\chi`` share `omega_c`, one ``K^{+}`` pass serves both the ``\phi`` and the ``\chi`` channel, so all three integrals are obtained in a single sweep over the kernel tables.

## Chemical potential
In the vDOS approximations the transition to the superconducting state shifts the chemical potential.
With `mu_flag = 1` (the default), IsoME determines ``\mu`` in every iteration by requiring that the electron number in the superconducting state matches the one in the normal state.
The root of ``N_e^{\text{nsc}}(\mu) - N_e^{\text{sc}}(\mu)`` is located with a Regula-Falsi method, starting from a ``\pm 50`` meV window around the current Fermi level that is shifted up or down in 50 meV steps until it brackets a sign change.
Only the ``Z``- and the ``\chi``-branch of the ``\varepsilon``-integral are needed for this, so the update is considerably cheaper than a full iteration: the ``W`` matrix product is never evaluated, although ``\phi_c`` still enters through the Lorentzian centres.

The ``\mu``-update is the part of the real-axis solver that reacts most sensitively to the grids and cutoffs; the [FAQ](@ref) lists the typical failure modes and how to recognize them.

## Self-consistency and convergence
The mixing is linear with an iteration-dependent factor — it starts at ``1``, decreases by ``0.05`` per iteration and settles at ``0.6`` from the ninth iteration on — applied to the primary variables ``Z``, ``\chi``, ``\phi_{ph}`` and ``\phi_c``.
A fixed factor can be enforced through `mixing_beta`.
The Coulomb contribution is ramped up over the first `nItFullCoul` iterations, which stabilizes the early iterations.

vDOS calculations are never started from scratch: the first temperature is seeded with a cDOS solution at ``\min(\mathrm{T}_C^{ML}/2,~T_1)``, where ``\mathrm{T}_C^{ML}`` is the machine-learning estimate printed with the Allen-Dynes table and ``T_1`` the first temperature solved (never below 0.5 K). Seeding below the target temperature is deliberate: the seed then carries a larger gap, which decays to zero if the temperature lies above ``\mathrm{T}_C``, whereas a seed with an almost vanishing gap may not grow past `minGap` in time below it. This cDOS run is attempted once; if it does not converge, the vDOS calculation starts from the BCS gap instead. Every subsequent temperature starts from the last converged solution.
``\mathrm{T}_C``-search mode uses the same bisection-and-fit strategy as the imaginary-axis solver, so ``\mathrm{T}_C`` ends up bracketed by the highest temperature whose gap edge stays above `minGap` and the lowest one without a solution.
Temperatures are stepped in whole kelvin, so the result is a 1 K bracket ``[T_{sc}, T_{nsc}]`` rather than a single number.

## Relevant input parameters
| Parameter | Role |
|:----------|:-----|
| `omega_c`, `domega` | Cutoff and step of the ``Z(\omega)``, ``\Delta(\omega)``, ``\chi(\omega)`` grid, of the ``\omega'``-tail, and of the kernel tables ``A(x)``, ``B(x)`` |
| `n_cheb` | Chebyshev points per pole in the head region of the ``\omega'``-grid |
| `depsilon` | Step of the ``\varepsilon``-grid (vDOS) |
| `encut` | Range of the ``\varepsilon``-integration |
| `mu_flag` | Whether the chemical potential is updated (vDOS only) |
| `mixing_beta` | Fixes the linear mixing factor instead of the default schedule |
| `nItFullCoul` | Iterations over which the Coulomb contribution is ramped up |
| `conv_thr`, `min_it` | Convergence threshold and minimum number of iterations |
| `minGap` | Gap below which a temperature counts as normal conducting |

See [Input](@ref) for defaults and types, and the [FAQ](@ref) for the practical consequences of these choices.

```@docs
RealAxisSolver
```
