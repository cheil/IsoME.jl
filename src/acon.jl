"""
    File containing Analytic continuation workflow:
        - Calculation of Green's functions
        - Nevanlinna
        - quasiparticle density of states

    Julia Packages:
        - Nevanlinna

    Comments:

"""


function acon(wsi, deltai, znormi, shifti, idx_ef)
    println("Analytic Continuation started")
    g11i, gauxi = calcGF(wsi, deltai[idx_ef, :], znormi, shifti)
    if lneva == 1
        ws, g11, gaux = NAC(wsi, real_c, g11i, gauxi)
        g12 = anomalousGF(g11, gaux)
        delta = Delta_from_GF(ws, g11, g12)
        dos_qp = qdos(ws, delta)#
    end

    lpade = 1
    if lpade == 1
        ws_pade, g11_pade, gaux_pade = Pade(wsi, nsiw, real_c, g11i, gauxi)
        g12_pade = anomalousGF(g11_pade, gaux_pade)
        delta_pade = Delta_from_GF(ws_pade, g11_pade, g12_pade)
        dos_qp_pade = qdos(ws_pade, delta_pade)
    end

    flag_A_figure = 1
    if flag_A_figure == 1
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
            A11 = -imag.(g11) / pi
            xlim_max = round(maximum(ws), RoundUp)
            xtick_val = 0:5:xlim_max
            ylim_max = round(maximum(A11), RoundUp)
            #ylim_max = round(maximum(A11) / 10 * 1.01, RoundUp) * 10
            #ylim_max = 5

            plot(ws, A11)
            xlims!(-xlim_max, xlim_max)
            ylims!(0, ylim_max)
            title!(material)
            xlabel!(L"\omega ~ \mathrm{(meV)}")
            ylabel!(L"\mathrm{A}^{\textrm{11}}(\omega) \;\;\; [\mathrm{eV^{-1}}]")
            savefig(outdir * "/A11_neva_T" * string(itemp) * ".pdf")

            A12 = -imag.(g12) / pi
            ylim_min = round(minimum(A12), RoundUp)
            ylim_max = round(maximum(A12), RoundUp)

            plot(ws, A12)
            xlims!(-xlim_max, xlim_max)
            ylims!(ylim_min, ylim_max)
            title!(material)
            xlabel!(L"\omega ~ \mathrm{(meV)}")
            ylabel!(L"\mathrm{A}^{\textrm{an}}(\omega) \;\;\; [\mathrm{eV^{-1}}]")
            savefig(outdir * "/A12_neva_T" * string(itemp) * ".pdf")
        end

        # Pade
        if lpade == 1
            A11_pade = -imag.(g11_pade) / pi   # this should be -imag.(g11)/pi, there is some sign issue...
            xlim_max = round(maximum(ws_pade), RoundUp)
            xtick_val = 0:10:xlim_max
            ylim_min = round(minimum(A11_pade), RoundUp)
            ylim_max = round(maximum(A11_pade), RoundUp)
            #ylim_max = round(maximum(A11) / 10 * 1.01, RoundUp) * 10
            #ylim_max = 5

            plot(ws_pade, A11_pade)
            xlims!(-xlim_max, xlim_max)
            ylims!(ylim_min, ylim_max)
            title!(material)
            xlabel!(L"\omega ~ \mathrm{(meV)}")
            ylabel!(L"\mathrm{A}^{\textrm{11}}(\omega) \;\;\; [\mathrm{eV^{-1}}]")
            savefig(outdir * "/A11_pade_T" * string(itemp) * ".pdf")

            A12_pade = -imag.(g12_pade) / pi   # this should be -imag.(g11)/pi, there is some sign issue...
            xlim_max = round(maximum(ws_pade), RoundUp)
            xtick_val = 0:10:xlim_max
            ylim_min = round(minimum(A12_pade), RoundUp)
            ylim_max = round(maximum(A12_pade), RoundUp)

            plot(ws_pade, A12_pade)
            xlims!(-xlim_max, xlim_max)
            ylims!(ylim_min, ylim_max)
            title!(material)
            xlabel!(L"\omega ~ \mathrm{(meV)}")
            ylabel!(L"\mathrm{A}^{\textrm{an}}(\omega) \;\;\; [\mathrm{eV^{-1}}]")
            savefig(outdir * "/A12_pade_T" * string(itemp) * ".pdf")
        end
    end

    if flag_qdos_figure == 1
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
            xlim_max = round(maximum(ws), RoundUp)
            xtick_val = 0:10:xlim_max
            #ylim_max = round(maximum(dos_qp) / 10 * 1.01, RoundUp) * 10
            ylim_max = round(maximum(dos_qp), RoundUp)

            plot(ws, dos_qp)
            xlims!(0, xlim_max)
            ylims!(0, ylim_max)
            title!(material)
            xlabel!(L"\omega ~ \mathrm{(meV)}")
            ylabel!(L"\mathrm{N}_{\mathrm{S}}(\omega)/\mathrm{N}_{\mathrm{F}}")
            savefig(outdir * "/qdos_neva_T" * string(itemp) * ".pdf")
        end

        # Pade
        if lpade == 1
            xlim_max = round(maximum(ws_pade), RoundUp)
            xtick_val = 0:10:xlim_max
            #ylim_max = round(maximum(dos_qp_pade) / 10 * 1.01, RoundUp) * 10
            ylim_max = round(maximum(dos_qp_pade), RoundUp)

            plot(ws_pade, dos_qp_pade)
            xlims!(0, xlim_max)
            ylims!(0, ylim_max)
            title!(material)
            xlabel!(L"\omega ~ \mathrm{(meV)}")
            ylabel!(L"\mathrm{N}_{\mathrm{S}}(\omega)/\mathrm{N}_{\mathrm{F}}")
            savefig(outdir * "/qdos_pade_T" * string(itemp) * ".pdf")
        end
    end
end

### Calculate Green's functions ###
function calcGF(wsi, deltai, znormi, shifti)
    """
    Calculate the normal and the auxiliary Green's function in imaginary space

    -------------------------------------------------------------------
    Input:
    wsi:        matsubara frequencies
    N:          number of matsubara frequencies
    deltai:     superconducting gap function in Matsubarsa frequency
    znormi:     renormalization function in Matsubara frequency
    shifti:     energy shift function in Matsubara frequency

    --------------------------------------------------------------------
    Output:
    g11i:       normal Green's function (11 component of Nambu)
                in imaginary space (complex)
    gauxi:      auxiliary Green's function in imaginary space (complex)

    --------------------------------------------------------------------
    """

    # println(size(wsi))
    # println(size(deltai))
    # println(size(znormi))
    # println(size(shifti))

    thetai = (wsi .* znormi).^2 .+ (shifti).^2 .+ (deltai .* znormi).^2
    g11i = -(im .* wsi  .* znormi .+ shifti) ./ thetai
    gauxi = -2 .* (im .* wsi  .* znormi .+ deltai .* znormi) ./ thetai

    return g11i, gauxi
end

### ACon ###
function Pade(wsi, nsiw, real_c, g11i, gauxi)
    """
    Calculate analytic continuation using Pade approximants

    --------------------------------------------------------------------
        Input:
        N:     number of matsubara frequencies
        g11i:     normal Green's function (11 component of Nambu) (complex)
        gauxi:    auxiliary Green's function                      (complex)

    --------------------------------------------------------------------
        Local:
        a11:        Pade coefficients
        aaux:       Pade coefficients

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
    if nsiw > 100
        N=100
    else
        N=nsiw
    end
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
    eta = 0.0005        # broadening parameter (w + im*eta)
    N_real = 1000       # dimension of array of output
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


function NAC(wsi, real_c, g11i, gauxi)
    """
    Use Nevanlinna package to analytically continue

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
    # needs "using Nevanlinna" at some point

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

### Nambu-Gor'kov components / Anomalous GF ###
function anomalousGF(g11, gaux)
    """
    Calculate anomalous Green's function

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

### Delta function from GF ###
function Delta_from_GF(ws, g11, g12)
    """
    Calculate Delta function from Green's functions

    --------------------------------------------------------------------
    Input:
    g11:        normal Green's function (11 component of Nambu)
                in real space (complex)
    g12:        anomalous Green's function (12 component of Nambu)
                in real space (complex)

    --------------------------------------------------------------------
    Output:
    delta:      superconducting gap function in real frequency
    """
    delta = (2 .* ws .* g12) ./ (g11 .- conj(reverse(g11)))
    return delta
end

### qdos ###
function qdos(ws, delta)
    """
    Calculate quasiparticle density of states

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

# plots for spectral function (Im of GFs) and qdos would be nice

