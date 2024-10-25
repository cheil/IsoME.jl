# IsoME

[![Build Status](https://github.com/cheil/IsoME.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/cheil/IsoME.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/cheil/IsoME.jl?svg=true)](https://ci.appveyor.com/project/cheil/IsoME-jl)


This short Julia codes solves the isotropic Migdal-Eliashberg equations, either within the constant DOS approximation or in the full-bandwidth (variable DOS) implementation. 

In both cases a file containing the Eliashberg spectral function alpha2F has to be provided.
For the variable DOS calculation a file with the electronic DOS and the respective value of the Fermi level has to be provided as well.

All input parameters and flags are set using the struct arguments.

The code is run by simply typing (with a working Julia setup in place):
julia main.jl




