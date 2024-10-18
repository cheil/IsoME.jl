"""
Solve the eliashberg equations for a given temperature. Cases:
    - 1) variable Dos with Weep
    - 2) variable Dos with Anderson pseudopotential
    - 3) constant Dos with Weep
    - 4) constatn Dos with Anderson pseudopotential

Julia Packages:
    - 

Comments:
    - 

"""

 
### Electron-Phonon Coupling ###
function calcLambda(itemp, M, a2f_omega, a2f)
    """
    Calculate electron-phonon coupling strength on imaginary axis
    based on ...
    λ(iω_n - iω_m) is calcualted for all pairs
    --> λ[1] = m-n = 0
        λ[2] = abs(m-n) = 1
        and so on

    -------------------------------------------------------------------
    Input:
        itemp:      temperature
        M:          number of highest matsubara frequency 
                    (attention: there is a 0'th frequency 
                    -> M=10 means there are 11 frequencies)
        a2f_omega:  energy in meV
        a2f_fine:   Eliashberg spectral function (1/meV) 

    --------------------------------------------------------------------
    Output:
        λ:          electron-phonon coupling for the first 2*M+2 
                    matsubara frequencies 
    
    --------------------------------------------------------------------
    Comments:
        - Calculate directly λ_n & λ_p?

    --------------------------------------------------------------------
    """    

    λ = zeros(2 * M + 2)
    for n in 1:2*M+2
        # n-th matsubara frequency
        omega = 2 * (n - 1) * π * kb * itemp

        # kernel 
        integrand = 2 * a2f .* a2f_omega ./ (omega .^ 2 .+ a2f_omega .^ 2)
        λ[n] = trapz(a2f_omega, integrand)
    end
    return λ
end


##### Eliashberg equations - Variable Dos - Weep #####
"""
    eliashberg_eqn(itemp, nsiw, wsi, ind_mat_freq, sparse_sampling_flag, lambdai, dosef, ndos, dos_en, dos, Weep, 
                   znormip, phiphip, phicip, shiftip, wgCoulomb, muintr)

Evaluate the isotropic Eliashberg equations, vDOS & Weep

"""
function eliashberg_eqn(itemp::Number, nsiw::Int64, wsi::Vector{Float64}, ind_mat_freq::Vector{Int64},
                        sparse_sampling_flag::Int64, lambdai::Vector{Float64}, dosef::Float64, ndos::Int64, 
                        dos_en::Vector{Float64}, dos::Vector{Float64}, Weep::Matrix{Float64}, znormip::Vector{Float64}, 
                        phiphip::Vector{Float64}, phicip::Vector{Float64}, shiftip::Vector{Float64}, wgCoulomb::Float64, 
                        muintr::Float64, idxEncut::Vector{Int64})
    """
    Evaluate the Isotropic Eliashberg equations for
        - variable Dos
        - Weep

    -------------------------------------------------------------------
    Input:
        itemp:          Temperature
        nsiw:           number of matsubara frequencies
        wsi:            matsubara frequencies
        ind_mat_freq:   indices of considered mat freq (sparse sampling)
        lambdai:        strength of el-ph interaction
        dosef:          Dos at fermi energy
        ndos:           number of energy points
        dos_en:         energies
        dos:            DoS
        znormip:        Z of previous iteration
        phiphip:        phi phonon of previous iteration
        phicip:         phi coloumb of previous iteration
        shiftip:        shift of previous iteration 
        wgCoulomb:         weighting of Weep 
        muintr:         updated chemical potential [optional]

    --------------------------------------------------------------------
    Output:
        data:   Vector containing the updated quantities
                    -
                    -
                    -

    --------------------------------------------------------------------
    Comments:
        - 
    --------------------------------------------------------------------
    """


    znormi  = zeros(nsiw)
    phici   = zeros(ndos)
    phiphi  = zeros(nsiw)
    shifti = zeros(nsiw)
    ckernel  = zeros(ndos,ndos)
    nsiw_vec = [1:nsiw;]

    # 1/theta matrix, NxM, i.e. ndos x nsiw
    theta_inv = 1 ./ ((znormip' .* wsi') .^ 2 .+ (phiphip' .+ phicip) .^ 2 .+ (dos_en .- muintr .+ shiftip') .^ 2)

    # matsubara sum, Nx1
    ckernel = sum(theta_inv .* (phiphip' .+ phicip), dims=2) - sum(phicip ./ (wsi' .^ 2 .+ (dos_en .- muintr) .^ 2 .+ phicip .^ 2), dims=2)

    # add analytic part, NxN
    ckernel = dos' .* wgCoulomb .* Weep .* (2 * kb * itemp .* ckernel' .+ 1 / 2 .* phicip' .* tanh.(1 / (2 * kb * itemp) .* sqrt.((dos_en' .- muintr)  .^ 2 .+ phicip' .^ 2)) ./ (sqrt.((dos_en' .- muintr) .^ 2 .+ phicip' .^ 2)))

    # integrate phic, Nx1    
    phici = -trapz(dos_en, ckernel)

    # z kernel, NxM
    zkernel =  dos .* theta_inv .* wsi' .* znormip'
    # ph kernel, NxM
    phkernel = dos .* theta_inv .* (phiphip' .+ phicip)
    # shift kernel, NxM
    shkernel = dos[idxEncut[1]:idxEncut[2]] .* theta_inv[idxEncut[1]:idxEncut[2],:] .* (dos_en[idxEncut[1]:idxEncut[2]] .- muintr .+ shiftip')

    # integrate over energy, gives vector of matsubara summands, Mx1
    ziwp = trapz(dos_en, zkernel')
    phiwp = trapz(dos_en, phkernel')
    shiwp = trapz(dos_en[idxEncut[1]:idxEncut[2]], shkernel')

    for iw in ind_mat_freq # loop over omega
        # Eq. (4.4) in Picket, PRB 26, 1186 (1982)
        tmp1 = lambdai[abs.(iw .- nsiw_vec).+1] 
        tmp2 = lambdai[iw.+nsiw_vec]
        lambdam = tmp1 .- tmp2
        lambdap = tmp1 .+ tmp2

        # Eqs. (4.1-4.3) in Picket, PRB 26, 1186 (1982) for FBW
        # using lambdam and lambdap the sum over |wp| < wscut 
        # is rewritten as a sum over iwp = 1, nsiw(itemp)
        znormi[iw] = znormi[iw] + dot(ziwp, lambdam)
        phiphi[iw] = phiphi[iw] + dot(phiwp, lambdap)
        shifti[iw] = shifti[iw] + dot(shiwp, lambdap)
    end 

    ### Interpolate remaining mat freg if sparse sampling ###
    if sparse_sampling_flag == 1
        znormi_sparse = filter(!iszero, znormi)
        phiphi_sparse = filter(!iszero, phiphi)
        shifti_sparse = filter(!iszero, shifti)

        znormi_itp = linear_interpolation(wsi[ind_mat_freq], znormi_sparse, extrapolation_bc=Line())
        phiphi_itp = linear_interpolation(wsi[ind_mat_freq], phiphi_sparse, extrapolation_bc=Line())
        shifti_itp = linear_interpolation(wsi[ind_mat_freq], shifti_sparse, extrapolation_bc=Line())

        znormi = znormi_itp[wsi[1:nsiw]]
        phiphi = phiphi_itp[wsi[1:nsiw]]
        shifti = shifti_itp[wsi[1:nsiw]]
    end

    # Eqs.(29)-(31) in Lee, Computational Materials (2023), 9:156
    znormi = 1.0 .+ kb * itemp .* znormi .* inv.(wsi) / dosef
    shifti = -kb * itemp .* shifti / dosef
    phiphi = kb * itemp .* phiphi / dosef

    data = Vector{Vector{Float64}}([znormi, phiphi, phici, shifti])

    return data
end


##### Eliashberg equations - Constant Dos - Weep #####
"""
    eliashberg_eqn( itemp, nsiw, wsi, ind_mat_freq, sparse_sampling_flag, lambdai, ndos, dos_en, 
                    dos, Weep, znormip, phiphip, phicip, idx_ef, wgCoulomb)

Evaluate the isotropic Eliashberg equations, cDOS & Weep
"""
function eliashberg_eqn(itemp::Number, nsiw::Int64, wsi::Vector{Float64}, ind_mat_freq::Vector{Int64}, 
                        sparse_sampling_flag::Int64, lambdai::Vector{Float64}, ndos::Int64, 
                        dos_en::Vector{Float64}, dos::Vector{Float64}, Weep::Matrix{Float64}, znormip::Vector{Float64}, 
                        phiphip::Vector{Float64}, phicip::Vector{Float64}, idx_ef::Int64, wgCoulomb::Float64)
    """
    Evaluate the Isotropic Eliashberg equations for
        - Constant Dos
        - Weep

    Following Eqs. (14)-(18) in Pellegrini, Phys- Mater. 5 024007(2022) 

    -------------------------------------------------------------------
    Input:
        itemp:          Temperature
        nsiw:           number of matsubara frequencies
        wsi:            matsubara frequencies
        ind_mat_freq:   indices of considered mat freq (sparse sampling)
        lambdai:        strength of el-ph interaction
        ndos:           number of energy points
        dos_en:         energies
        dos:            DoS
        znormip:        Z of previous iteration
        phiphip:        phi phonon of previous iteration
        phicip:         phi coloumb of previous iteration
        wgCoulomb:         weighting of Weep 

    --------------------------------------------------------------------
    Output:
        data:   Vector containing the updated quantities
                    -
                    -
                    -

    --------------------------------------------------------------------
    Comments:
        - 
    --------------------------------------------------------------------
    """

    # Solve Eliashberg equations for current iteration
    znormi  = zeros(nsiw)
    phici   = zeros(ndos)
    phiphi  = zeros(nsiw)
    ckernel  = zeros(ndos,ndos)

    nsiw_vec = [1:nsiw;]

    # denominator z & ph, Mx1
    normzph = 1 ./ sqrt.((znormip .* wsi) .^ 2 .+ (phiphip .+ phicip[idx_ef]) .^ 2) #+wz²+phi⁰ #[ef]
    
    # 1/theta, NxM
    theta_inv = 1 ./ ((znormip' .* wsi') .^ 2 .+ (phiphip' .+ phicip) .^ 2 .+ (dos_en .^ 2)) #theta[iwn](e)vector in epsilon

    # vector of matsubara summands, Mx1
    ziwp = normzph .* wsi .* znormip
    phiwp = normzph .* (phiphip .+ phicip[idx_ef])

    # perform matsubara sum, Nx1
    ckernel = sum(theta_inv .* (phiphip' .+ phicip) .- (phicip ./ (wsi' .^ 2 .+ dos_en .^ 2 .+ phicip .^ 2)), dims=2)

    # add analytic part and multiply with Weep, NxN
    ckernel = dos' .* wgCoulomb .* Weep .* (2 * kb * itemp .* ckernel' .+ 1 / 2 .* phicip' .* tanh.(1 / (2 * kb * itemp) .* sqrt.(dos_en' .^ 2 .+ phicip' .^ 2)) ./ (sqrt.(dos_en' .^ 2 .+ phicip' .^ 2)))

    # integrate phi_c
    phici = -trapz(dos_en, ckernel)

    for iw in ind_mat_freq # loop over omega
        # Eq. (4.4) in Picket, PRB 26, 1186 (1982)
        tmp1 = lambdai[abs.(iw .- nsiw_vec).+1] 
        tmp2 = lambdai[iw.+nsiw_vec]
        lambdam = tmp1 .- tmp2
        lambdap = tmp1 .+ tmp2

        znormi[iw] = znormi[iw] + dot(ziwp, lambdam)
        phiphi[iw] = phiphi[iw] + dot(phiwp, lambdap)
    end # iw = 1:M+1

    # sparse sampling 
    if sparse_sampling_flag == 1
        znormi_sparse = filter(!iszero, znormi)
        phiphi_sparse = filter(!iszero, phiphi)

        znormi_itp = interpolate((wsi[ind_mat_freq],), znormi_sparse, Gridded(Linear()))
        phiphi_itp = interpolate((wsi[ind_mat_freq],), phiphi_sparse, Gridded(Linear()))

        znormi = znormi_itp[wsi[1:nsiw]]
        phiphi = phiphi_itp[wsi[1:nsiw]]
    end 

    # Eqs. (14)-(18) in Pellegrini, Phys- Mater. 5 024007(2022) 
    znormi = 1.0 .+ π * kb * itemp .* znormi ./ wsi
    phiphi = π * kb * itemp .* phiphi

    data = Vector{Vector{Float64}}([znormi, phiphi, phici])

    return data
end


##### Eliashberg equations - Variable Dos - mu* #####
"""
    eliashberg_eqn( itemp, nsiw, wsi, ind_mat_freq, sparse_sampling_flag, lambdai, dos_en, 
                    dos, dosef, znormip, deltaip, shiftip, muc, muintr, wgCoulomb)

Evaluate the isotropic Eliashberg equations, variable dos, mu*
"""
function eliashberg_eqn(itemp::Number, nsiw::Int64, wsi::Vector{Float64}, ind_mat_freq::Vector{Int64}, 
                        sparse_sampling_flag::Int64, lambdai::Vector{Float64}, dos_en::Vector{Float64}, 
                        dos::Vector{Float64}, dosef::Float64, znormip::Vector{Float64}, deltaip::Vector{Float64}, 
                        shiftip::Vector{Float64}, muc::Float64, muintr::Float64, wgCoulomb::Float64, idxEncut::Vector{Int64})
    """
    Evaluate the Isotropic Eliashberg equations for
        - Variable Dos
        - no Weep
    
    Following Lee, npj Computational Materials volume 9, 156 (2023)
        - Eqs.(29) - (31)

    -------------------------------------------------------------------
    Input:
        itemp:          Temperature
        nsiw:           number of matsubara frequencies
        wsi:            matsubara frequencies
        ind_mat_freq:   indices of considered mat freq (sparse sampling)
        lambdai:        strength of el-ph interaction
        dosef:          Dos at fermi energy
        dos_en:         energies
        dos:            DoS
        znormip:        Z of previous iteration
        deltaip:        delta of previous iteration
        shiftip:        shift of previous iteration 
        muc_ME:         mu star for Migdal-Eliashberg
        muintr:         updated chemical potential [optional]

    --------------------------------------------------------------------
    Output:
        data:   Vector containing the updated quantities
                    - znormi
                    - deltai
                    - shifti

    --------------------------------------------------------------------
    Comments:
        - (ϵ, iωₙ) <--> (dos_en, wsi')
    --------------------------------------------------------------------
    """

    deltai = zeros(nsiw)
    znormi = zeros(nsiw)
    shifti = zeros(nsiw)
    nsiw_vec = [1:nsiw;]

    theta_inv = 1 ./ (znormip' .^ 2 .* (wsi' .^ 2 .+ deltaip' .^ 2) .+ (dos_en .- muintr .+ shiftip') .^ 2)
    zkernel = dos .* theta_inv .* wsi' .* znormip'
    dekernel = dos .* theta_inv .* deltaip' .* znormip'
    shkernel = dos[idxEncut[1]:idxEncut[2]] .* theta_inv[idxEncut[1]:idxEncut[2],:] .* (dos_en[idxEncut[1]:idxEncut[2]] .- muintr .+ shiftip')

    ziwp = trapz(dos_en, zkernel') 
    deiwp = trapz(dos_en, dekernel') 
    shiwp = trapz(dos_en[idxEncut[1]:idxEncut[2]], shkernel') 


    for iw in ind_mat_freq # loop over omega
        tmp1 = lambdai[abs.(iw .- nsiw_vec).+1] 
        tmp2 = lambdai[iw.+nsiw_vec]
        lambdam = tmp1 .- tmp2
        lambdap = tmp1 .+ tmp2
        # Eqs. (4.1-4.3) in Picket, PRB 26, 1186 (1982) for FBW
        # using lambdam and lambdap the sum over |wp| < wscut 
        # is rewritten as a sum over iwp = 1, nsiw(itemp)
        znormi[iw] = znormi[iw] + dot(ziwp, lambdam)
        deltai[iw] = deltai[iw] + dot(deiwp, lambdap .- 2 * wgCoulomb * muc)
        shifti[iw] = shifti[iw] + dot(shiwp, lambdap)
    end # iw = 1:M+1

    # gives different result? (e.g. 1 K difference in Tc for LaBeH8, especially shift is different)
    if sparse_sampling_flag == 1
        znormi_sparse = filter(!iszero, znormi)
        deltai_sparse = filter(!iszero, deltai)
        shifti_sparse = filter(!iszero, shifti)

        znormi_itp = linear_interpolation(wsi[ind_mat_freq], znormi_sparse, extrapolation_bc=Line())
        deltai_itp = linear_interpolation(wsi[ind_mat_freq], deltai_sparse, extrapolation_bc=Line())
        shifti_itp = linear_interpolation(wsi[ind_mat_freq], shifti_sparse, extrapolation_bc=Line())

        znormi = znormi_itp[wsi[1:nsiw]]
        deltai = deltai_itp[wsi[1:nsiw]]
        shifti = shifti_itp[wsi[1:nsiw]]
    end

    # Eqs.(29)-(31) in Lee, npj Computational Materials volume 9, 156 (2023)
    znormi = 1.0 .+ kb * itemp .* znormi .* inv.(wsi) / dosef
    deltai = kb * itemp .* deltai .* inv.(znormi) / dosef
    shifti = -kb * itemp .* shifti / dosef

    #plot(wsi, shifti)
    #savefig("shifti_temp.png")

    data = Vector{Vector{Float64}}([znormi, deltai, shifti])

    return data
end


##### Eliashberg equations - Constant Dos - mu* #####
function eliashberg_eqn(itemp::Number, nsiw::Int64, wsi::Vector{Float64}, ind_mat_freq::Vector{Int64},
                        sparse_sampling_flag::Int64, lambdai::Vector{Float64}, deltaip::Vector{Float64}, 
                        muc::Float64, wgCoulomb)
    """
    Evaluate the Isotropic Eliashberg equations for
        - Constant Dos
        - no Weep

    Following Margine & Giustino, Phys. Rev. B 87, 024505
        - Eq. (34) - (35)

    -------------------------------------------------------------------
    Input:
        itemp:          Temperature
        nsiw:           number of matsubara frequencies
        wsi:            matsubara frequencies
        ind_mat_freq:   indices of considered mat freq (sparse sampling)
        lambdai:        strength of el-ph interaction
        dos_en:         energies
        znormip:        Z of previous iteration
        phiphip:        phi phonon of previous iteration
        muc_ME:         mu star for Migdal-Eliashberg

    --------------------------------------------------------------------
    Output:
        data:   Vector containing the updated quantities

    --------------------------------------------------------------------
    Comments:
        - 
    --------------------------------------------------------------------
    """


    deltai = zeros(nsiw)
    znormi = zeros(nsiw)
    nsiw_vec = [1:nsiw;]

    normzph = 1 ./ sqrt.(wsi .^ 2 .+ deltaip .^ 2) #wz^2+phi0

    ziwp = normzph .* wsi 
    deiwp = normzph .* deltaip

    for iw in ind_mat_freq # loop over omega
        tmp1 = lambdai[abs.(iw .- nsiw_vec).+1] 
        tmp2 = lambdai[iw.+nsiw_vec]
        lambdam = tmp1 .- tmp2
        lambdap = tmp1 .+ tmp2
        # Eqs. (4.1-4.3) in Picket, PRB 26, 1186 (1982) for FBW
        # using lambdam and lambdap the sum over |wp| < wscut 
        # is rewritten as a sum over iwp = 1, nsiw(itemp)
        znormi[iw] = znormi[iw] + dot(ziwp, lambdam)
        deltai[iw] = deltai[iw] + dot(deiwp, lambdap .- 2 * wgCoulomb * muc)
    end # iw = 1:M+1


    #gives different result? (e.g. 1 K difference in Tc for LaBeH8, especially shift is different)
    if sparse_sampling_flag == 1
        znormi_sparse = filter(!iszero, znormi)
        deltai_sparse = filter(!iszero, deltai)

        znormi_itp = linear_interpolation(wsi[ind_mat_freq], znormi_sparse, extrapolation_bc=Line())
        deltai_itp = linear_interpolation(wsi[ind_mat_freq], deltai_sparse, extrapolation_bc=Line())

        znormi = znormi_itp[wsi[1:nsiw]]
        deltai = deltai_itp[wsi[1:nsiw]]
    end

    # Eqs.(34)-(35) in Margine and Giustino, PRB 87, 024505 (2013)
    znormi = 1.0 .+ π * kb * itemp .* znormi .* inv.(wsi)
    deltai = π * kb * itemp .* deltai .* inv.(znormi)


    data = Vector{Vector{Float64}}([znormi, deltai])

    return data
end



function initSparseSampling(beta, omega_c, M)
    """
    Initialize sparse matsubara basis (IR-Basis)

    -------------------------------------------------------------------
    Input:
        beta:       inverse temperature
        omega_c:    omega cutoff
        M:          number mat freq

    --------------------------------------------------------------------
    Output:
        ind_mat_freq:   indices of relevant matsubara frequencies
    
    --------------------------------------------------------------------
    Comments:
        - 
    --------------------------------------------------------------------
    """


    IRbasis = FiniteTempBasis(Fermionic(), beta, omega_c)
    ir_mat  = MatsubaraSampling(IRbasis; positive_only=true)
    ir_mat  = SparseIR.value.(ir_mat.ωn, beta)
    ir_indices_str = MatsubaraSampling(IRbasis; positive_only=true)
    ir_indices     = zeros(Int, size(ir_mat))

    for i_ir ∈ 1:length(ir_mat)
        dummy = ir_indices_str.sampling_points[i_ir]
        ir_indices[i_ir] = (dummy.n + 1) / 2 #+1
    end

    ind_mat_freq = ir_indices[ir_indices.<M+2] #+2

    if ind_mat_freq[end] != M+1 # last mat freq needed for interpolation; 
        ind_mat_freq = [ind_mat_freq; M+1] 
    end

    return ind_mat_freq

end


