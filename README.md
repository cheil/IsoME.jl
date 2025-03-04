# IsoME

[![Build Status](https://github.com/cheil/IsoME.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/cheil/IsoME.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/cheil/IsoME.jl?svg=true)](https://ci.appveyor.com/project/cheil/IsoME-jl)


This Julia codes solves the isotropic Migdal-Eliashberg equations, either within the constant DOS approximation, the full-bandwidth (variable DOS) implementation or with the full static coulomb interaction $W(\epsilon, \epsilon')$.

In all cases a file containing the Eliashberg spectral function alpha2F has to be provided.
For the variable DOS calculation the electronic DOS and the respective value of the Fermi level is needed as well.
In the most general case an additional file containing the screened Coulomb interaction W is required.

All input parameters and flags are set using the custom struct `arguments`.

For a more thorough documentation please refer to [Docs](https://cheil.github.io/IsoME.jl/).





