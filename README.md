# IsoME

[![Build Status](https://github.com/cheil/IsoME.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/cheil/IsoME.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/cheil/IsoME.jl?svg=true)](https://ci.appveyor.com/project/cheil/IsoME-jl)


This short Julia codes solves the isotropic Migdal-Eliashberg equations, either within the constant DOS approximation or in the full-bandwidth (variable DOS) implementation. 

In both cases a file containing the Eliashberg spectral function alpha2F has to be provided.
For the variable DOS calculation a file with the electronic DOS and the respective value of the Fermi level has to be provided as well.

All input parameters and flags are set using the struct arguments.

The code is run by simply typing (with a working Julia setup in place):
julia main.jl


# ToDo's
- convergence test for mu*_ME 
- Write documentation / create package
    * When creating a package a certain structure has to be obeyed
    * The documentation can be included to the package (.gitlab-ci.yml)
         + use "Documenter" packager 
           --> minimal example https://gitlab.com/gitlab-examples/julia
           --> documentation https://documenter.juliadocs.org/stable/man/guide/#Package-Guide
    * Adapt docu of functions https://docs.julialang.org/en/v1/manual/documentation/
    * add .toml file to documentation
- convergence threshold??
- Improve Mixing
- Error when starting from terminal??
- Tc search mode: For a given material the code should be able to find the Tc w/o running through all temperatures, the user may give an initial guess for the Tc
    * Manual mode: user specifies all temperatures 
    * Auto mode:   User specifies just one initial guess for T
- termination criterion: if Delta0 < 0.1 for e.g. first 20 iterations, then force a higher initial guess for gap0??
- Gap increases sometimes when other dos is used
    * Change Tc search mode s.t. if gap increases manual mode is used
    * Occurs if fsthick is too small (~ < 30 eV)


--- Discarded / On hold ---
- Shift energies in vDos to ef=0?
    * Problem: in each iteration epsilon - mu is calculated --> epsilon + ef - mu should work
    * Did not work because problem with mu update --> why?
- replace trapz by simps? --> Only if there is an official julia package 
- Broyden mixing:
    * Initial values  (cDos W)                
        + dZ = zeros(nBroyd, nsiw)
        + dZ[end,:] = znormi
        + dZ_F = zeros(nBroyd, nsiw)
        + G_Z = I*0.5
        + dph = zeros(nBroyd, nsiw)
        + dph[end,:] = phiphi
        + dph_F = zeros(nBroyd, nsiw)
        + G_ph = I*0.5
        + dc = zeros(nBroyd, ndos)
        + dc[end,:] = phici
        + dc_F = zeros(nBroyd, ndos)
        + G_c = I*0.5 
    * mixing  (cDos W) 
        + znormi, dZ, dZ_F = broyden_mixing(i_it, znormip, dZ, new_data[1]-znormip, dZ_F, I, I)
        + phiphi, dph, dph_F = broyden_mixing(i_it, phiphip, dph, new_data[2]-phiphip, dph_F, G_ph, I, I)
        + phici, dc, dc_F = broyden_mixing(i_it, phicip, dc, new_data[3]-phicip, dc_F, G_c, I, I)



---- Done ----
- Error handle, error file
- readin: cut off zero of DOS at beginning and end?
- Format log-file analog to console output?
    * write elapsed time to console
- specify units as string instead of numerical value, e.g. "eV" --> create flag from that | Auto extraction of flag from header of each file
- Let user specify spin factor in Dos
- CLean up flag zoo
    * Optional: Input flags
    * Discarded: integration, broyden, 
- IR basis only for T < 1; no user input, internal
- fsthick no user parameters
- write formulas into paper
- Speed up no Weep code + rewrite in terms of Delta (Following Margine & Giustino)
- Regular falsi instead of bisection?
- use sparse also without Weep?
- single function for weep and no weep (e.g. merge solve_eliashberg_noWeep and solve_eliashberg)
- include fsthick also for Weep? 
- Find good sampling points for interpolation of ph_ph/Z !!
- check number of needed epsilons for meaningful integration --> INTERPOLATE !!
- phi_c times tanh (-> error in paper?)
- how to initialize phi_c? (phi_c(ef) can not be zero, otherwise we divide by zero for phi_c)
- Interpolate outliers in Weep (see plot_W.jl)
- Interpolate kernel and use more integration points !!
- Flag include W / mu* 
- mu update (does not converge sometimes) --> Bisection method
- improve perfomance (one for instead of two over iw and iwp) 
- check speed-up with IR basis (bottelneck: calculation of phi_c takes 10-100 times longer -> 10x speedup due to IR basis is not significant; e.g. one iteration takes 15 sec without IR and 14,9 with IR)

*1000 to read in Weep (done)
check if phi_ph is used instead of phi? (done?)
Delta_0_new - correct to take phi[1,1]?
eq. 16 paper: include density of states at N(e')
check dimension mismatch? -> debug
2* phi c! compare paper Picket, PRB 26, 1186 (1982) and Sanna
check/rewrite input done by hand (index of fermi energy in weep, energy window for epsion values)
plot W 3d - symmetric

# FAQ
- Convergence parameter:
    * omega_c
    * fsthick
- Material converges to completely unexpected Tc
    * Check if all unit flags are set correctly / have been extracted correctly 
- The code uses the wrong units
    * check if a wrong unit occurs in the header of any of the input files 
- Gap increases with Temperature
    * DOS may not converged 
    * fsthick is maybe too low / DOS is given in a window too narrow --> this is in particular a problem when a other dos is used 
- mu converges to strange values
    * Try to increase omega_c / decrease fsthick
- Specifications input files:
    * Weep: data has to be in third column
    * Dos: energies must be in first, dos in second column
    * a2f: first column energies, 2:end a2f for different smearings
- Tc search mode: unexpected T < 0 K
    * check if gap increases with T
- Specified units are not recognized
    * Units are case sensitive
    * Currently available are: eV, meV, THz, Ry

