# IsoME

[![Build Status](https://github.com/cheil/IsoME.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/cheil/IsoME.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://cheil.github.io/IsoME.jl/)

This Julia codes solves the isotropic Migdal-Eliashberg equations, either within the constant DOS approximation, the full-bandwidth (variable DOS) implementation or with the full static coulomb interaction $W(\epsilon, \epsilon')$.

In all cases a file containing the Eliashberg spectral function $\alpha^2F$ has to be provided.
For the variable DOS calculation the electronic DOS $N(\varepsilon)$ and the respective value of the Fermi level is needed as well.
In the most general case an additional file containing the screened Coulomb interaction $W(\varepsilon,\varepsilon')$ is required.

All input parameters and flags are set using the custom struct `arguments`.

## Installation

Julia 1.10 or higher is required. IsoME is a registered package:

```julia
julia> using Pkg
julia> Pkg.add("IsoME")
```

## The two solvers

Both solvers take the same `arguments` instance and the same input files.

```julia
using IsoME

inp = arguments(a2f_file = "Nb.a2f", outdir = "Nb_run")

EliashbergSolver(inp)   # imaginary (Matsubara) axis - fast, robust, the choice for Tc searches
RealAxisSolver(inp)     # direct real-frequency axis - spectral quantities without analytic continuation
```

Output is written to `outdir`, which defaults to `IsoME/` in the current working directory.
An existing directory is never written into; a run counter is appended instead (`IsoME_1/`, `IsoME_2/`, ...).

Runnable examples for all approximations are in [`test/Nb/examples.jl`](test/Nb/examples.jl).

## Documentation

For a more thorough documentation please refer to [Docs](https://cheil.github.io/IsoME.jl/).
Breaking changes between releases are listed in [`CHANGELOG.md`](CHANGELOG.md).

## Citing

If you use IsoME in your work, please cite

> E. Kogler, D. Spath, R. Lucrezi, H. Mori, Z. Zhu, Z. Li, E. R. Margine and C. Heil,
> *IsoME: Streamlining high-precision Eliashberg calculations*,
> Computer Physics Communications **315**, 109720 (2025).
> [doi:10.1016/j.cpc.2025.109720](https://doi.org/10.1016/j.cpc.2025.109720)

A BibTeX entry is in [`CITATION.bib`](CITATION.bib).
