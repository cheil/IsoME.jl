"""
    file containing the real axis eliashberg equations
"""

"""
    realEliashbergEq(beta, znormip, phi_ph_ip, phi_c_ip, shiftip, ws, w_static, dosef,
                     fermi_level, wgCoulomb, electronic_spec)

Real axis Eliashberg equations in the vDOS+W approximation.
The ε-grid, the DOS and W(ε,ε′) are not passed separately: they enter only through the
piecewise-linear moments M0/M1 and M0W/M1W carried in `electronic_spec`.
"""
function realEliashbergEq(beta::Float64, znormip::Vector{ComplexF64}, phi_ph_ip::Vector{ComplexF64}, phi_c_ip::Vector{ComplexF64}, shiftip::Vector{ComplexF64},
                          ws::WPrimeWorkspace, w_static::AbstractVector,
                          dosef::Float64, fermi_level::Float64, wgCoulomb::Float64, electronic_spec::Tuple)

    w_prime = ws.wp_full

    # interpolate Z,χ,ϕ onto ω'-integration grid
    Z_itp = linear_interpolation(w_static, znormip, extrapolation_bc=Flat())
    shift_itp = linear_interpolation(w_static, shiftip, extrapolation_bc=Flat())

    Z_ongrid = Z_itp.(w_prime)
    shift_ongrid = shift_itp.(w_prime) .- fermi_level
    # φ(ω,ε) = φ_ph(ω) + φ_c(ε): only φ_ph needs the ω→ω' interpolation; φ_c lives on the ε-grid.
    phi_ph_ongrid = linear_interpolation(w_static, phi_ph_ip, extrapolation_bc=Flat()).(w_prime)

    # ------------- ε-integration ------------- #
    integrands, coulomb_spectral = epsilon_helpers_vDOSW(electronic_spec, Z_ongrid, phi_ph_ongrid, phi_c_ip, shift_ongrid, w_prime)

    # Z-integrand must be positive (causality), flip the others accordingly
    # Note: Minus sign because split of -Θ to (ε+χ +- ε_p) (eq. (51))
    FLIP = ifelse.(integrands[1] .<= 0, 1, -1)
    g_z = -abs.(integrands[1])
    g_phi = -(FLIP .* integrands[2])
    g_chi = -(FLIP .* integrands[3])

    # ------------- Ω & ω'-integration ------------- #
    # O(N)-memory kernel route, one pass: K⁻ g_z plus K⁺ g_phi and K⁺ g_chi, all on w_static.
    Iz, Iphi, Ichi = kernel_omega_integral(ws, w_static, g_z, g_phi, g_chi)

    # Coulomb term: integrate N(ε')W(ε,ε') with the same piecewise-linear spectral
    # quadrature (coulomb_spectral is ndos × M). The dosef prefactor cancels.
    wgt_full = vcat(ws.wgt_head, ws.wgt_tail)
    coulomb_weight = wgt_full .* FLIP .* tanh.(beta .* w_prime ./ 2)
    phi_c_val = (wgCoulomb / π) .* (coulomb_spectral * coulomb_weight)   # φ_c(ε), length ndos

    Zval = 1 .+ Iz ./ (w_static .* (π * dosef))
    phi_ph_val = Iphi ./ (π * dosef)                                     # φ_ph(ω), length nω
    shift_val = -Ichi ./ (π * dosef)

    # φ(ω,ε) = φ_ph(ω) + φ_c(ε) is kept in separable form; no nω×ndos matrix is built.
    # Δ is derived from φ and Z in the solver after mixing; hand back the two φ vectors here.
    return Zval, shift_val, phi_ph_val, phi_c_val
end


"""
    realEliashbergEq(mu_star, beta, znormip, phiip, shiftip, ws, w_static, dosef, wgCoulomb,
                     fermi_level, electronic_spec)

Real axis Eliashberg equations in the vDOS+μ approximation. 
"""
function realEliashbergEq(mu_star::Float64, beta::Float64, znormip::Vector{ComplexF64}, phiip::Vector{ComplexF64},
                          shiftip::Vector{ComplexF64}, ws::WPrimeWorkspace, w_static::AbstractVector,
                          dosef::Float64, wgCoulomb::Float64,
                          fermi_level::Float64, electronic_spec::Tuple)

    w_prime = ws.wp_full

    # interpolate Z,χ,ϕ onto ω'-integration grid (φ is handed in directly)
    Z_itp = linear_interpolation(w_static, znormip, extrapolation_bc=Flat())
    phi_itp = linear_interpolation(w_static, phiip, extrapolation_bc=Flat())
    shift_itp = linear_interpolation(w_static, shiftip, extrapolation_bc=Flat())

    Z_ongrid = Z_itp.(w_prime)
    phi_ongrid = phi_itp.(w_prime)
    shift_ongrid = shift_itp.(w_prime) .- fermi_level

    # ------------- ε-integration ------------- #
    integrands = epsilon_helpers(electronic_spec, Z_ongrid, phi_ongrid, shift_ongrid, w_prime)

    # Z-integrand must be positive (causality), flip the others accordingly
    # Note: Minus sign because split of -Θ to (ε+χ +- ε_p) (eq. (51))
    FLIP = ifelse.(integrands[1] .<= 0, 1, -1)
    g_z = -abs.(integrands[1])
    g_phi = -(FLIP .* integrands[2])
    g_chi = -(FLIP .* integrands[3])

    # ------------- ω'-integration ------------- #
    # O(N)-memory kernel route, one pass: K⁻ g_z plus K⁺ g_phi and K⁺ g_chi, all on w_static.
    Iz, Iphi, Ichi = kernel_omega_integral(ws, w_static, g_z, g_phi, g_chi)

    # μ* Coulomb term
    coulomb = wprime_trapz(ws, g_phi .* tanh.(beta .* w_prime ./ 2))

    Zval = 1 .+ Iz ./ (w_static .* (π * dosef))
    phi_val = (Iphi .- wgCoulomb .* mu_star .* coulomb) ./ (π * dosef)
    shift_val = -Ichi ./ (π * dosef)

    return Zval, shift_val, phi_val
end




"""
     realEliashbergEq(mu_star, beta, deltaip, ws, w_static)

Real axis Eliashberg equations in cDOS+μ approximation
"""
function realEliashbergEq(mu_star::Float64, beta::Float64, deltaip::Vector{ComplexF64},
                          ws::WPrimeWorkspace, w_static::AbstractVector, wgCoulomb::Float64)

    Delta_func = linear_interpolation(w_static, deltaip, extrapolation_bc=Flat())

    w_prime = ws.wp_full
    Delta_ongrid = Delta_func.(w_prime)
    sqrt_eval = @. sqrt(w_prime^2 - Delta_ongrid^2)

    g_z = @. real(w_prime / sqrt_eval)          # Z   integrand density (∝ quasiparticle DOS)
    g_phi = @. real(Delta_ongrid / sqrt_eval)   # Δ·Z integrand density

    # O(N)-memory kernel route: materialized K head (gemv) + Toeplitz/Hankel tail (A + f·B on the fly)
    Iz, Iphi = kernel_omega_integral(ws, w_static, g_z, g_phi)

    # μ* Coulomb term (scalar, same for all ω)
    coulomb = wprime_trapz(ws, g_phi .* tanh.(beta .* w_prime ./ 2))

    Zval = 1 .- Iz ./ w_static
    Delta = (Iphi .- wgCoulomb .* mu_star .* coulomb) ./ Zval

    return Zval, Delta
end



###################################################
# ------------------- Helpers ------------------- #
###################################################
"""
    epsilon_helpers(electronic_spec, Z_ongrid, phi_ongrid, shift_ongrid, w_prime; idx_skip = nothing)

Sovle the epsilon integral in the vDOS+μ approximation
Returns the ω'-integrands
"""
function epsilon_helpers(electronic_spec, Z_ongrid, phi_ongrid, shift_ongrid, w_prime; idx_skip = nothing)
    eps_1, eps_2, dos_1, _, dε, ddos, _, _ = electronic_spec

    # assuming a linear form of the dos; ddos/dε are the precomputed differences
    M1 = transpose(ddos ./ dε)
    M0 = transpose(dos_1) .- transpose(eps_1) .* M1
    
    # Helper functions for the ε-integration
    ε_p = sqrt.(w_prime.^2 .*Z_ongrid.^2 .- phi_ongrid.^2)       
    Rplus = real.(shift_ongrid .+ ε_p)
    Rminus = real.(shift_ongrid .- ε_p)
    Iplus = imag.(shift_ongrid .+ ε_p)
    Iminus = imag.(shift_ongrid .- ε_p)

    # M0- and M1-weighted ε-sums of the Lorentzian moments (pp and mm branches)
    lorentzian_moments_p = lorentzian_moments_weighted(eps_1, eps_2, Rplus, Iplus, M0, M1)
    lorentzian_moments_m = lorentzian_moments_weighted(eps_1, eps_2, Rminus, Iminus, M0, M1)

    # Evaluate ε-integrals to get ω'-integrands 
    integrands = Vector{Vector{Float64}}(undef, 3)
    for (idx, g) in enumerate((w_prime .* Z_ongrid, phi_ongrid, shift_ongrid))
        if ~isnothing(idx_skip) && any(idx .== idx_skip)
            continue
        end
        
        integrands[idx] = eval_spectral_integrals(g, lorentzian_moments_p, lorentzian_moments_m, ε_p, idx == 3)
    end

    return integrands
end


"""
    spectral_interval_moments(εl, εr, A, B)

The six combinations of the Lorentzian moments Jⁿ = ∫ εⁿ/((ε+A)²+B²) dε over [εl, εr]
that the spectral coefficients actually need, for one (ω′, ε-interval) pair:

    Tₙ = Jⁿ + A·Jⁿ⁻¹   (n = 1,2,3)      Uₙ = B·Jⁿ   (n = 0,1,2)

Jⁿ is never returned on its own, and that is the whole point. Splitting the partial
fraction into an odd and an even part about ε = -A,

    Im[ (g/2ε_p) / (ε+χ±ε_p) ] = I_g·(ε+R±)/D  -  I±·R_g/D ,

shows that the arctangent enters only through those two groups: the odd part contributes
T = Jⁿ + A·Jⁿ⁻¹, in which the arctangent cancels identically (T₁ = H/2, pure log), and the
even part contributes U = B·Jⁿ, in which the explicit damping factor cancels the 1/B of the
collapsing Lorentzian. So the assembled integrand is smooth as B → 0 even when the pole
ε = -A sits inside the interval - the Sokhotski-Plemelj limit, whose δ-function pickup is
finite.

Forming J⁰ ~ π/B first and multiplying by B afterwards, as the code used to, throws that
away: the intermediate overflows the mantissa long before B reaches zero, and at B = 0 the
product is 0·∞ = NaN. Here the bounded primitive

    W = B·J⁰ = atan2(B·Δε, B² + P),   P = (A+εr)(A+εl),   |W| ≤ π

is built directly, via atan(u) - atan(v) = atan2(u-v, 1+uv) rescaled by B² > 0 (quadrant
preserving, so an identity). No division by B occurs anywhere. The same regrouping fixes a
second, far more common failure: for the vast majority of intervals the pole lies outside
(P > 0) and Jⁿ is perfectly finite, but the old endpoint subtraction cancelled two O(1/B)
numbers to produce it and lost every digit. 
See `Overleaf RealMe Derivations/lorentzian_moments.tex`.
"""
@inline function spectral_interval_moments(εl, εr, A, B)
    A2 = A^2
    B2 = B^2

    Δε = εr - εl
    xl = A + εl
    xr = A + εr
    # < 0 exactly when the pole ε = -A lies inside [εl, εr]
    P = xr * xl

    W = atan(B * Δε, B2 + P)                    # B·J⁰, bounded by π, exact at B = 0
    H = log((xr^2 + B2) / (xl^2 + B2))          # one log of the ratio, not a difference

    T1 = 0.5 * H                                # J¹ + A·J⁰ - the arctangent drops out
    T2 = Δε - 0.5 * A * H - B * W               # J² + A·J¹
    T3 = 2 * A * B * W + 0.5 * (A2 - B2) * H + 0.5 * Δε * (εr + εl - 2 * A)

    U0 = W                                      # B·J⁰
    U1 = 0.5 * B * H - A * W                    # B·J¹
    U2 = B * Δε - A * B * H + (A2 - B2) * W     # B·J²

    return T1, T2, T3, U0, U1, U2
end

"""
    spectral_C0C1(Rg, Ig, Ta_p, Tb_p, Ua_p, Ub_p, Ta_m, Tb_m, Ua_m, Ub_m)

Scalar version of `eval_spectral_integral_coefficients` for one (ω′, ε-interval) pair,
in terms of the grouped moments of `spectral_interval_moments`:

    C0 = I_g·(Tₐ⁺ - Tₐ⁻) - R_g·(Uₐ⁺ - Uₐ⁻),    C1 = I_g·(T_b⁺ - T_b⁻) - R_g·(U_b⁺ - U_b⁻)

The slope (g_slope) contribution has the same functional form with the moments shifted by
one order, so it is obtained by passing (T2, T3, U1, U2) instead of (T1, T2, U0, U1).
R±/I± no longer appear here - they are already folded into T and U.
"""
@inline function spectral_C0C1(Rg, Ig, Ta_p, Tb_p, Ua_p, Ub_p, Ta_m, Tb_m, Ua_m, Ub_m)
    C0 = Ig * (Ta_p - Ta_m) - Rg * (Ua_p - Ua_m)
    C1 = Ig * (Tb_p - Tb_m) - Rg * (Ub_p - Ub_m)
    return C0, C1
end

"""
    epsilon_helpers_vDOSW(electronic_spec, Z_ongrid, phi_ph_ongrid, phi_c, shift_ongrid, w_prime)
    epsilon_causality_vDOSW(electronic_spec, Z_ongrid, phi_ph_ongrid, phi_c, shift_ongrid, w_prime)

Fused route for the vDOS+W spectral integrals: instead of materialising the
Lorentzian-moment matrices (and S, P, R±, I± as full ω′×ε matrices),
loop over the ε-intervals and accumulate the M0/M1-weighted ε-sum directly
into length-M vectors for the Z, φ and χ integrands.
Only the φ-branch coefficients C0/C1 are stored as matrices, because the Coulomb
part couples them to M0W/M1W over the *output* ε-grid via a GEMM.

`epsilon_causality_vDOSW` is the μ-update's variant: the same sweep without the φ branch and
without the Coulomb matrices, since neither is needed to conserve charge, and it runs on every
root-finder evaluation of `diff_Ne_realAxis`. Both dispatch into `_epsilon_helpers_vDOSW` with
the φ flag as a `Val`, so the branch is a compile-time constant and each specialisation is the
machine code of a hand-written loop - measured neutral to within 1% against the two separate
functions this replaced, on identical data and bit-identical output.

W(ε,ε′) is not an argument: it enters only through M0W/M1W in `electronic_spec`.
"""
epsilon_helpers_vDOSW(electronic_spec, Z_ongrid, phi_ph_ongrid, phi_c, shift_ongrid, w_prime) =
    _epsilon_helpers_vDOSW(electronic_spec, Z_ongrid, phi_ph_ongrid, phi_c, shift_ongrid, w_prime, Val(true))

epsilon_causality_vDOSW(electronic_spec, Z_ongrid, phi_ph_ongrid, phi_c, shift_ongrid, w_prime) =
    _epsilon_helpers_vDOSW(electronic_spec, Z_ongrid, phi_ph_ongrid, phi_c, shift_ongrid, w_prime, Val(false))

function _epsilon_helpers_vDOSW(electronic_spec, Z_ongrid, phi_ph_ongrid, phi_c, shift_ongrid, w_prime,
                                ::Val{WITH_PHI}) where {WITH_PHI}
    eps_1, eps_2, dos_1, _, dε, ddos, M0W, M1W = electronic_spec

    M = length(w_prime)
    Nint = length(dε)

    z_int = zeros(Float64, M)
    phi_int = zeros(Float64, M)
    chi_int = zeros(Float64, M)

    # φ-branch coefficients: needed as full matrices for the Coulomb GEMM below
    C0phi = Matrix{Float64}(undef, WITH_PHI ? M : 0, WITH_PHI ? Nint : 0)
    C1phi = Matrix{Float64}(undef, WITH_PHI ? M : 0, WITH_PHI ? Nint : 0)

    wZ = w_prime .* Z_ongrid

    @inbounds for jε in 1:Nint
        xl = eps_1[jε]
        xr = eps_2[jε]
        invdε = 1.0 / dε[jε]
        m1 = ddos[jε] * invdε
        m0 = dos_1[jε] - xl * m1

        # φ(ω,ε) = φ_ph(ω) + φ_c(ε): the ε-slope Φ1 = φ_c'(ε) is ω-independent
        Φ1 = (phi_c[jε+1] - phi_c[jε]) * invdε
        inv_scale = 1 / (1 + Φ1^2)

        for iw in 1:M
            Φ0 = phi_ph_ongrid[iw] + phi_c[jε] - xl * Φ1
            sh = shift_ongrid[iw]
            wz = wZ[iw]

            S = (sh + Φ0 * Φ1) * inv_scale
            P = sqrt((wz^2 - sh^2 - Φ0^2) * inv_scale + S^2)

            Sp = S + P
            Sm = S - P
            Rp = real(Sp); Ip = imag(Sp)
            Rm = real(Sm); Im_ = imag(Sm)

            T1p, T2p, T3p, U0p, U1p, U2p = spectral_interval_moments(xl, xr, Rp, Ip)
            T1m, T2m, T3m, U0m, U1m, U2m = spectral_interval_moments(xl, xr, Rm, Im_)

            c = inv_scale / (2 * P)

            # Z branch: g = ωZ/scale, no slope
            gz = wz * c
            C0, C1 = spectral_C0C1(real(gz), imag(gz),
                                   T1p, T2p, U0p, U1p, T1m, T2m, U0m, U1m)
            z_int[iw] += m0 * C0 + m1 * C1

            # φ branch: g = Φ0/scale, slope = Φ1/scale (μ-update does not need it)
            if WITH_PHI
                gphi = Φ0 * c
                gs = Φ1 * c
                C0, C1 = spectral_C0C1(real(gphi), imag(gphi),
                                       T1p, T2p, U0p, U1p, T1m, T2m, U0m, U1m)
                D0, D1 = spectral_C0C1(real(gs), imag(gs),
                                       T2p, T3p, U1p, U2p, T2m, T3m, U1m, U2m)
                C0 += D0
                C1 += D1
                C0phi[iw, jε] = C0
                C1phi[iw, jε] = C1
                phi_int[iw] += m0 * C0 + m1 * C1
            end

            # χ branch: g = χ/scale, slope = 1/scale
            gchi = sh * c
            C0, C1 = spectral_C0C1(real(gchi), imag(gchi),
                                   T1p, T2p, U0p, U1p, T1m, T2m, U0m, U1m)
            D0, D1 = spectral_C0C1(real(c), imag(c),
                                   T2p, T3p, U1p, U2p, T2m, T3m, U1m, U2m)
            chi_int[iw] += m0 * (C0 + D0) + m1 * (C1 + D1)
        end
    end

    if WITH_PHI
        coulomb_spectral = Matrix{Float64}(undef, size(M0W, 1), M)
        mul!(coulomb_spectral, M0W, transpose(C0phi))
        mul!(coulomb_spectral, M1W, transpose(C1phi), 1.0, 1.0)
        return [z_int, phi_int, chi_int], coulomb_spectral
    else
        return z_int, chi_int
    end
end



"""
    lorentzian_moments_weighted(x0, x1, A, B, M0, M1)

Evaluate the grouped Lorentzian moments of `spectral_interval_moments` together with the
weighting factors M0/M1, summing directly over ε.
A, B are R+/I+ or R-/I-. B is ω′-dependent only, so it factors out of the ε-sum and the
grouping survives it.

Only 8 of the 12 possible weight/moment pairs are ever read downstream, so only those are
accumulated: `eval_spectral_integrals` keeps C0 from the M0-weighted call, which uses
(T1, T2, U0, U1), and C1 from the M1-weighted call, which uses (T2, T3, U1, U2). The other
four (M0·T3, M0·U2, M1·T1, M1·U0) were computed and discarded.
"""
function lorentzian_moments_weighted(x0, x1, A::AbstractVector, B::AbstractVector, M0, M1)
    xl = vec(x0); xr = vec(x1); m0 = vec(M0); m1 = vec(M1)
    Nε = length(xl)
    M  = length(A)
    M0_T1 = zeros(M); M0_T2 = zeros(M); M0_U0 = zeros(M); M0_U1 = zeros(M)
    M1_T2 = zeros(M); M1_T3 = zeros(M); M1_U1 = zeros(M); M1_U2 = zeros(M)
    @inbounds for jε in 1:Nε
        l = xl[jε]; r = xr[jε]; w0 = m0[jε]; w1 = m1[jε]
        Δε = r - l
        @simd for iw in 1:M
            a = A[iw]; b = B[iw]
            b2 = b^2; a2 = a^2
            # grouped endpoint differences - see spectral_interval_moments
            pl = a + l; pr = a + r
            W = atan(b * Δε, b2 + pr * pl)
            H = log((pr^2 + b2) / (pl^2 + b2))

            T1 = 0.5 * H
            T2 = Δε - 0.5 * a * H - b * W
            T3 = 2 * a * b * W + 0.5 * (a2 - b2) * H + 0.5 * Δε * (r + l - 2 * a)
            U0 = W
            U1 = 0.5 * b * H - a * W
            U2 = b * Δε - a * b * H + (a2 - b2) * W

            M0_T1[iw] += w0 * T1; M0_T2[iw] += w0 * T2
            M0_U0[iw] += w0 * U0; M0_U1[iw] += w0 * U1
            M1_T2[iw] += w1 * T2; M1_T3[iw] += w1 * T3
            M1_U1[iw] += w1 * U1; M1_U2[iw] += w1 * U2
        end
    end

    return (M0_T1, M0_T2, M0_U0, M0_U1, M1_T2, M1_T3, M1_U1, M1_U2)
end


"""
    eval_spectral_integrals(g, lorentzian_moments_p, lorentzian_moments_m, ε_p, shift_int; g_slope)

Evaluate spectral integrals for wZ, delta, chi based on the precomputed lorentzian integrals
The lorentzian moments are already weighted by M0/M1 and the ε-sumamtion is done
"""
function eval_spectral_integrals(g, lorentzian_moments_p, lorentzian_moments_m, ε_p, shift_int::Bool=false; g_slope=nothing)
    slope = (shift_int && isnothing(g_slope)) ? one.(ε_p) : g_slope

    M0_T1_p, M0_T2_p, M0_U0_p, M0_U1_p, M1_T2_p, M1_T3_p, M1_U1_p, M1_U2_p = lorentzian_moments_p
    M0_T1_m, M0_T2_m, M0_U0_m, M0_U1_m, M1_T2_m, M1_T3_m, M1_U1_m, M1_U2_m = lorentzian_moments_m

    # C0 from the M0-weighted sums, C1 from the M1-weighted ones; the slope part is the same
    # functional form with the moments shifted by one order (T1,U0 -> T2,U1 -> T3,U2)
    C0s = eval_spectral_integral_coefficients(g, ε_p, M0_T1_p, M0_T2_p, M0_U0_p, M0_U1_p,
                                                      M0_T1_m, M0_T2_m, M0_U0_m, M0_U1_m; g_slope=slope)
    C1s = eval_spectral_integral_coefficients(g, ε_p, M1_T2_p, M1_T3_p, M1_U1_p, M1_U2_p,
                                                      M1_T2_m, M1_T3_m, M1_U1_m, M1_U2_m; g_slope=slope)
    return C0s .+ C1s
end


"""
    eval_spectral_integral_coefficients(g, ε_p, Ta_p, Tb_p, Ua_p, Ub_p, Ta_m, Tb_m, Ua_m, Ub_m; g_slope)

Vectorised counterpart of `spectral_C0C1`: one coefficient from one pair of grouped moments
per branch. The slope term of an affine numerator g(ε) = g + g_slope·ε uses the same form one
order up, which is why (Ta, Ua) and (Tb, Ub) are passed together. The vDOS+μ χ integral
corresponds to g_slope = 1.
"""
function eval_spectral_integral_coefficients(g, ε_p, Ta_p, Tb_p, Ua_p, Ub_p,
                                                     Ta_m, Tb_m, Ua_m, Ub_m; g_slope=nothing)
    Rg = real(g ./ (2 * ε_p))
    Ig = imag(g ./ (2 * ε_p))

    C = @. Ig * (Ta_p - Ta_m) - Rg * (Ua_p - Ua_m)

    if !isnothing(g_slope)
        Rg = real(g_slope ./ (2 * ε_p))
        Ig = imag(g_slope ./ (2 * ε_p))
        @. C += Ig * (Tb_p - Tb_m) - Rg * (Ub_p - Ub_m)
    end

    return C
end
