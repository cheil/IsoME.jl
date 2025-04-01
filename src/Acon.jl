"""
    File containing Analytic continuation workflow:
        - Calculation of Green's functions
        - Nevanlinna
        - quasiparticle density of states

    Julia Packages:
        - Nevanlinna

    Comments:

"""

# TODO:
#   - flags as input
#   - real_c as input
function acon(inp, itemp, wsi, nsiw, deltai, znormi, log_file; idx_ef=-1, shifti = 0)

    real_c = 100.0

    # sinnvoller machen
    if inp.include_Weep == 1
        deltai = deltai[idx_ef, :]
    end

    g11i, gauxi = calcGF(wsi, deltai, znormi, shifti)
    lneva = 0
    if lneva == 1
        ws, g11, gaux = NAC(wsi, real_c, g11i, gauxi)
        g12 = anomalousGF(g11, gaux)
        delta = Delta_from_GF(ws, g11, g12)
        dos_qp = qdos(ws, delta)#

        # spectral function
        A11 = -imag.(g11) / pi
        A12 = -imag.(g12) / pi

    end

    lpade = 1
    if lpade == 1
        #=
        ------- direct pade of GF, less stable -------
        ws_pade, g11_pade, gaux_pade = Pade(wsi, nsiw, real_c, g11i, gauxi)
        g12_pade = anomalousGF(g11_pade, gaux_pade)
        delta_pade = Delta_from_GF(ws_pade, g11_pade, g12_pade)
        dos_qp_pade = qdos(ws_pade, delta_pade)

        # spectral function
        A11_pade = -imag.(g11_pade) / pi
        A12_pade = -imag.(g12_pade) / pi
        ----------------------------------------------
        =#

        delta_pade, ws_pade = Pade_separate(wsi, nsiw, real_c, deltai)
    end

    ### save outputs ###
    folder = inp.outdir*"ACON/"
    if ~isdir(folder)
        mkdir(folder)
    end

    if lneva == 1
        saveACON(itemp, folder, ws, delta, A11, A12, "neva")
    end

    if lpade == 1
        # saveACON(itemp, folder, ws_pade, delta_pade, A11_pade, A12_pade, "pade") # direct Pade of GF
        saveACON(itemp, folder, ws_pade, delta_pade, "pade")
    end    

    # make plots
    flag_figureACON = 1
    if flag_figureACON == 1
        plot_font = "Computer Modern"
        default(
            fontfamily=plot_font,
            linewidth=2,
            framestyle=:box,
            label=nothing,
            grid=false
        )

        # NAC
        if lneva == 1
            plotSpectralFunction(itemp, ws, A11, A12, folder, inp.material)

            plotQDOS(itemp, ws, dos_qp, folder, inp.material)

            plotGap(itemp, ws, delta, folder, inp.material)
        end

        # Pade
        if lpade == 1
            #=
            ------- direct pade of GF, less stable -------
            plotSpectralFunction(itemp, ws_pade, A11_pade, A12_pade, folder, inp.material, "pade")

            plotQDOS(itemp, ws_pade, dos_qp_pade, folder, inp.material, "pade")
            ----------------------------------------------
            =#

            plotGap(itemp, ws_pade, delta_pade, folder, inp.material, "pade")
        end
    end

    printTee(log_file, "Analytic Continuation finished\n")
end



"""
    calcGF(wsi, deltai, znormi, shifti)

Calculate the normal and the auxiliary Green's function in imaginary space
"""
function calcGF(wsi, deltai, znormi, shifti)

    thetai = (wsi .* znormi).^2 .+ (shifti).^2 .+ (deltai .* znormi).^2
    g11i = -(im .* wsi  .* znormi .+ shifti) ./ thetai
    gauxi = -2 .* (im .* wsi  .* znormi .+ deltai .* znormi) ./ thetai

    return g11i, gauxi
end

"""
    Pade(wsi, nsiw, real_c, g11i, gauxi)

Calculate analytic continuation of Green's function using Pade approximants
"""
function Pade(wsi, nsiw, real_c, g11i, gauxi)
    """
    --------------------------------------------------------------------
        Input:
        g11i:     normal Green's function (11 component of Nambu) (complex)
        gauxi:    auxiliary Green's function                      (complex)

    --------------------------------------------------------------------
        Output:
        ws:         real frequency vector (not broadened)
        g11:        normal Green's function (11 component of Nambu)
                    in real space (complex)
        gaux:       auxiliary Green's function in real space (complex)

    --------------------------------------------------------------------
    """

    # pade coeff
    # Pade collapses with to many matsubara points:
    N = minimum([100, nsiw])

    f11i = Array{Complex}(undef, N, N)
    fauxi = Array{Complex}(undef, N, N)
    a11 = Array{Complex}(undef, N)
    aaux = Array{Complex}(undef, N)
    for p in 1:N
        if p == 1
            f11i[p, :] = g11i[1:N]
            fauxi[p, :] = gauxi[1:N]
        else
            for i in p:N
                tmp1 = f11i[p-1, p-1] / f11i[p-1, i]
                tmp2 = f11i[p-1, i] / f11i[p-1, i]
                f11i[p, i] = (tmp1 - tmp2) / (wsi[i]*im - wsi[p-1]*im)
                tmp1 = fauxi[p-1, p-1] / fauxi[p-1, i]
                tmp2 = fauxi[p-1, i] / fauxi[p-1, i]
                fauxi[p, i] = (tmp1 - tmp2) / (wsi[i]*im - wsi[p-1]*im)
            end
        end
        a11[p] = f11i[p, p]
        aaux[p] = fauxi[p, p]
        #println("a11=",a11[p])
    end
    # Check whether a[p] is NaN
    ar = real.(a11)
    ai = imag.(a11)
    if any(isnan.(ar)) || any(isnan.(ai))
        println("One or more Pade coefficients are NaN")
    end
    ar = real.(aaux)
    ai = imag.(aaux)
    if any(isnan.(ar)) || any(isnan.(ai))
        println("One or more Pade coefficients are NaN")
    end
    # pade eval
    eta = 0.001        # broadening parameter (w + im*eta)
    N_real = 5000       # dimension of array of output
    ws = LinRange(-real_c, real_c, N_real)
    
    A = Array{Complex}(undef, N_real, N+1)
    B = Array{Complex}(undef, N_real, N+1) # pade approximants
    A[:,1] .= 0
    A[:,2] .= a11[1]
    B[:,1] .= 1
    B[:,2] .= 1

    for i in 3:N+1
        # Version with broadening:
        #A[:,i] = A[:,i-1] .+ (ws.+im*eta .- wsi[i-2]*im) .* a11[i-1] .* A[:,i-2]
        #B[:,i] = B[:,i-1] .+ (ws.+im*eta .- wsi[i-2]*im) .* a11[i-1] .* B[:,i-2]
        A[:,i] = A[:,i-1] .+ (ws .- wsi[i-2]*im) .* a11[i-1] .* A[:,i-2]
        B[:,i] = B[:,i-1] .+ (ws .- wsi[i-2]*im) .* a11[i-1] .* B[:,i-2]
    end
    #println("final A=",A[1:10,N-50])
    #println("final B=",B[:,N]) #these are NaN!

    # Check whether a[p] is NaN
    Ar = real.(A)
    Ai = imag.(A)
    Br = real.(B)
    Bi = imag.(B)
    if any(isnan.(Ar)) || any(isnan.(Ai)) || any(isnan.(Br)) || any(isnan.(Bi))
        println("One or more Pade approximant is NaN")
        return
    end

    g11 = A[:,N+1]./B[:,N+1] # elementwise for vectors of ws, otherwise scalar

    A[:,1] .= 0
    A[:,2] .= aaux[1]
    B[:,1] .= 1
    B[:,2] .= 1
    for i in 3:N+1
        A[:,i] = A[:,i-1] .+ (ws .- wsi[i-2]*im) .* aaux[i-1] .* A[:,i-2]
        B[:,i] = B[:,i-1] .+ (ws .- wsi[i-2]*im) .* aaux[i-1] .* B[:,i-2]
    end
    gaux = A[:,N+1]./B[:,N+1]
    return ws, g11, gaux
end



"""
    Pade(wsi, nsiw, real_c, g11i, gauxi)

Calculate analytic continuation of complex function using Pade approximants
"""
function Pade_separate(wsi, nsiw, real_c, gi)
    """
    Calculate analytic continuation of a complex function Pade approximants

    --------------------------------------------------------------------
        Input:
        wsi:      Matsubara frequencies (real, meV)
        nsiw:     Number of Matsubara frequencies
        real_c:   Real frequency cutoff [meV]
        gi:       Complex function gi(wsi) to be analytically continued

    --------------------------------------------------------------------
        Local:
        a:        Pade coefficients
        A, B:     Pade approximants

    --------------------------------------------------------------------
        Output:
        ws:       Real frequency vector (not broadened)
        g:        Analytic continuation of complex function g(ws)

    --------------------------------------------------------------------
    """
    # pade coeff
    # Pade collapses with to many matsubara points:
    #println("Start Pade_seperate")
    if nsiw > 200
        N=200
    else
        N=nsiw
    end
    #N = nsiw
    f = Array{Complex}(undef, N, N)
    a = Array{Complex}(undef, N)
    #println("Arrays defined")
    for p in 1:N
        if p == 1
            f[p,1:N] = gi[1:N]
        else
            for i in p:N
                tmp1 = f[p-1, p-1] / f[p-1, i]
                tmp2 = f[p-1, i] / f[p-1, i]
                f[p, i] = (tmp1 - tmp2) / (wsi[i]*im - wsi[p-1]*im)
            end
        end
        a[p] = f[p, p]
        #println("a=",a[p])
    end
    #println("Finished a's")
    # Check whether a[p] is NaN
    ar = real.(a)
    ai = imag.(a)
    if any(isnan.(ar)) || any(isnan.(ai))
        println("One or more Pade coefficients are NaN")
    end
    # pade eval
    eta = 0.001        # broadening parameter (w + im*eta)
    N_real = 2000       # dimension of array of output
    ws = LinRange(-real_c, real_c, N_real)
    
    A = Array{Complex}(undef, N_real, N+1)
    B = Array{Complex}(undef, N_real, N+1) # pade approximants
    A[:,1] .= 0 + 0*im
    A[:,2] .= a[1]
    B[:,1] .= 1 + 0*im
    B[:,2] .= 1 + 0*im
    for i in 3:N+1
        # Version with broadening:
        A[:,i] = A[:,i-1] .+ ((ws.+im*eta) .- wsi[i-2]*im) .* a[i-1] .* A[:,i-2]
        B[:,i] = B[:,i-1] .+ ((ws.+im*eta) .- wsi[i-2]*im) .* a[i-1] .* B[:,i-2]
        # A[:,i] = A[:,i-1] .+ (ws .- wsi[i-2]*im) .* a[i-1] .* A[:,i-2]
        # B[:,i] = B[:,i-1] .+ (ws .- wsi[i-2]*im) .* a[i-1] .* B[:,i-2]
    end
    # println("final A=",A[1:10,N-50])
    # println("final B=",B[:,N])
    #println("Finished A,B")

    # Check whether A nad B is NaN
    Ar = real.(A)
    Ai = imag.(A)
    Br = real.(B)
    Bi = imag.(B)
    if any(isnan.(Ar)) || any(isnan.(Ai)) || any(isnan.(Br)) || any(isnan.(Bi))
        println("One or more Pade approximants is NaN")
    end

    g = A[:,N+1]./B[:,N+1] # elementwise for vectors of ws, otherwise scalar
    #println("Finished g")
    return g, ws
end


"""
    NAC(wsi, real_c, g11i, gauxi)

Use Nevanlinna package to analytically continue
"""
function NAC(wsi, real_c, g11i, gauxi)
    """
    --------------------------------------------------------------------
    Input:
    wsi:        Matsubara frequencies
    real_c:     energy cutoff of real axis [meV]
    g11i:       normal Green's function (11 component of Nambu)
                in imaginary space (complex)
    gauxi:      auxiliary Green's function in imaginary space (complex)

    --------------------------------------------------------------------
    Output:
    ws:         real frequency vector (not broadened)
    g11:        normal Green's function (11 component of Nambu)
                in real space (complex)
    gaux:       auxiliary Green's function in real space (complex)

    --------------------------------------------------------------------
    """

    # these are required, and could be user inputs:
    eta = 0.0005         # broadening parameter (w + im*eta)
    N_real = 4000       # dimension of array of output

    # these are fixed:
    norm = 1.0          # sum rule: Integral of Spectral function normalized to 1
    norm_aux = 2.0      # sum rule: Integral of Spectral function normalized to 2

    # these might as well not exist, if optimization=false is not done, but the function
    # still needs values for some reason:
    H_max = 50          # cutoff for Hardy basis
    lambda = 1e-4       # regularization parameter
    iter_tol = 500      # upper bound of iteration

    sol = Nevanlinna.NevanlinnaSolver(wsi*im, g11i, N_real, real_c, eta, norm, H_max, iter_tol, lambda, verbose=false, optimization=false)
    sol_aux = Nevanlinna.NevanlinnaSolver(wsi*im, gauxi, N_real, real_c, eta, norm_aux, H_max, iter_tol, lambda, verbose=false, optimization=false)

    ws = real.(sol.reals.freq) # has negative and positive values, imag part is broadening
    g11 = -sol.reals.val        # -conj(g11(-w))=g22(w)
    gaux = -sol_aux.reals.val   # NAC package for some reason flips sign
    return ws, g11, gaux
end

"""
    anomalousGF(g11, gaux)

Calculate anomalous Green's function
"""
function anomalousGF(g11, gaux)
    """

    --------------------------------------------------------------------
    Input:
    g11:        normal Green's function (11 component of Nambu)
                in real space (complex, must have negative real values) 
    gaux:       auxiliary Green's function in real space (complex)

    --------------------------------------------------------------------
    Output:
    g12:        anomalous Green's function (12 component of Nambu)
                in real space (complex)

    --------------------------------------------------------------------
    """
    g12 = 0.5 .* (gaux .- g11 .+ conj(reverse(g11)))
    return g12
end


"""
    Delta_from_GF(ws, g11, g12)

Calculate gap function from Green's functions
"""
function Delta_from_GF(ws, g11, g12)

    delta = (2 .* ws .* g12) ./ (g11 .- conj(reverse(g11)))

    return delta
end

### Add Z and chi as in Daniels paper eq. (27), (28)

"""
    qdos(ws, delta)

Calculate quasiparticle density of states
"""
function qdos(ws, delta)
    """
    --------------------------------------------------------------------
    Input:
    ws:         real frequency vector
    delta:      superconducting gap function in real frequency

    --------------------------------------------------------------------
    Output:
    dos_qp:     quasiparticle density of states in BCS-limit

    """

    eta = 0.0005
    omega = ws .+ im * eta

    dos_qp = real.(omega ./ sqrt.(omega.^2 .- delta.^2))

    return dos_qp
end


### Plots ###
"""

Plot the spectral function on the real axis
"""
function plotSpectralFunction(itemp, ws, A11, A12, folder, material, mode="neva")

    xlim_max = round(maximum(ws), RoundUp)
    xtick_val = 0:5:xlim_max
    ylim_max = round(maximum(A11), RoundUp)
    #ylim_max = round(maximum(A11) / 10 * 1.01, RoundUp) * 10
    #ylim_max = 5

    plot(ws, A11)
    xlims!(-xlim_max, xlim_max)
    ylims!(0, ylim_max)
    if material != "Material"
        title!(inp.material)
    end
    xlabel!(L"\omega ~ \mathrm{(meV)}")
    ylabel!(L"\mathrm{A}^{\textrm{11}}(\omega) \;\;\; [\mathrm{eV^{-1}}]")
    savefig(folder * "/A11_"*mode*"_T" * string(itemp) * ".pdf")

    ylim_min = round(minimum(A12), RoundUp)
    ylim_max = round(maximum(A12), RoundUp)

    plot(ws, A12)
    xlims!(-xlim_max, xlim_max)
    ylims!(ylim_min, ylim_max)
    if material != "Material"
        title!(inp.material)
    end
    xlabel!(L"\omega ~ \mathrm{(meV)}")
    ylabel!(L"\mathrm{A}^{\textrm{an}}(\omega) \;\;\; [\mathrm{eV^{-1}}]")
    savefig(folder * "/A12_"*mode*"_T" * string(itemp) * ".pdf") 

end

"""

plot the QDOS in the BCS limit on the real axis
"""
function plotQDOS(itemp, ws, dos_qp, folder, material, mode="neva")
    xlim_max = round(maximum(ws), RoundUp)
    xtick_val = 0:10:xlim_max
    #ylim_max = round(maximum(dos_qp) / 10 * 1.01, RoundUp) * 10
    ylim_max = round(maximum(dos_qp), RoundUp)

    plot(ws, dos_qp)
    xlims!(0, xlim_max)
    ylims!(0, ylim_max)
    if material != "Material"
        title!(inp.material)
    end
    xlabel!(L"\omega ~ \mathrm{(meV)}")
    ylabel!(L"\mathrm{N}_{\mathrm{S}}(\omega)/\mathrm{N}_{\mathrm{F}}")
    savefig(folder* "/qdos_"*mode*"_T" * string(itemp) * ".pdf")
end


"""

plot the superconducting gap function on the real axis
"""
function plotGap(itemp, ws, delta, folder, material, mode= "neva")

    xlim_max = round(maximum(ws), RoundUp)
    xtick_val = 0:10:xlim_max
    #ylim_max = round(maximum(delta), RoundUp)

    plot(ws, real(delta), label="real", color = :red)
    plot!(ws, imag(delta), label="imag", color = :blue)
    xlims!(0, xlim_max)
    #ylims!(0, ylim_max)
    if material != "Material"
        title!(inp.material)
    end
    xlabel!(L"\omega ~ \mathrm{(meV)}")
    ylabel!(L"\Delta(\omega) ~ \mathrm{(meV)}")
    savefig(folder* "/gap_"*mode*"_T" * string(itemp) * ".pdf")
end


"""
    saveSelfEnergyComponents(inp, iwn, Delta, Z; epsilon=nothing,  chi=nothing, phiph=nothing, phic=nothing)

Save the the real frequency quantities
"""
function saveACON(itemp, folder, ws, delta, A11, A12, mode="neva")

    # ACON
    open(folder*"ACON_"*mode*"_"*string(itemp)*"K.dat", "a") do io
        write(io, "#  ω / meV       A11(ω) / 1       A12(ω) / 1        Δ(ω) / meV\n")
        writedlm(io, zip(ws, A11, A12, delta), '\t')
    end

end


"""
    saveSelfEnergyComponents(inp, iwn, Delta, Z; epsilon=nothing,  chi=nothing, phiph=nothing, phic=nothing)

Save the the real frequency gap
"""
function saveACON(itemp, folder, ws, delta, mode="neva")

    # ACON
    open(folder*"ACON_"*mode*"_"*string(itemp)*"K.dat", "a") do io
        write(io, "#  ω / meV          Δ(ω) / meV\n")
        writedlm(io, zip(ws, delta), '\t')
    end

end