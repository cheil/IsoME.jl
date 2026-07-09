"""
    file containing the real axis eliashberg equations
"""


"""
    realEliashbergEq(beta, znormip, phiphip, shiftip, ws, w_static, w_static_chi, dosef,
                     epsilon, dos, Weep, idx_ef, fermi_level, wgCoulomb, electronic_spec)

Real axis Eliashberg equations in the vDOS+W approximation
"""
function realEliashbergEq(beta::Float64, znormip::Vector{ComplexF64}, phiphip::Matrix{ComplexF64}, shiftip::Vector{ComplexF64},
                          ws::WPrimeWorkspace, w_static::AbstractVector,
                          w_static_chi, dosef::Float64, epsilon::Vector{Float64}, dos::Vector{Float64}, Weep::Matrix{Float64},
                          idx_ef::Int64, fermi_level::Float64, wgCoulomb::Number, electronic_spec::Tuple)

    if size(Weep, 1) != length(epsilon) || size(Weep, 2) != length(epsilon)
        error("The W(ε,ε′) matrix must be defined on the same energy grid as the DOS for real-axis vDOS+W calculations.")
    end

    w_prime = ws.wp_full

    # interpolate Z,χ,ϕ onto ω'-integration grid
    Z_itp = linear_interpolation(w_static, znormip, extrapolation_bc=Flat())
    shift_itp = linear_interpolation(w_static_chi, shiftip, extrapolation_bc=Flat())

    Z_ongrid = Z_itp.(w_prime)
    shift_ongrid = shift_itp.(w_prime) .- fermi_level
    phi_ongrid = interpolate_phi_matrix(w_static, phiphip, w_prime)

    # ------------- ε-integration ------------- #
    integrands, coulomb_spectral = eval_spectral_and_coulomb_vDOS_W(electronic_spec, epsilon, dos, Weep, Z_ongrid, phi_ongrid, shift_ongrid, w_prime)

    # Z-integrand must be positive (causality), flip the others accordingly
    FLIP = ifelse.(integrands[1] .<= 0, 1, -1)
    g_z = -abs.(integrands[1])
    g_phi = -(FLIP .* integrands[2])
    g_chi = -(FLIP .* integrands[3])

    # ------------- Ω & ω'-integration ------------- #
    Iz = kernel_omega_integral(ws.Km_head, ws.Km_tail, ws, g_z)
    Iphi = kernel_omega_integral(ws.Kp_head, ws.Kp_tail, ws, g_phi)
    Ichi = kernel_omega_integral(ws.Kpchi_head, ws.Kpchi_tail, ws, g_chi)

    # Coulomb term: integrate N(ε')W(ε,ε') with the same piecewise-linear spectral
    # quadrature (coulomb_spectral is ndos × M). The dosef prefactor cancels.
    wgt_full = vcat(ws.wgt_head, ws.wgt_tail)
    coulomb_weight = wgt_full .* FLIP .* tanh.(beta .* w_prime ./ 2)
    phi_c_val = (wgCoulomb / π) .* (coulomb_spectral * coulomb_weight)   # length ndos

    Zval = 1 .+ Iz ./ (w_static .* (π * dosef))
    phi_ph_val = Iphi ./ (π * dosef)
    shift_val = -Ichi ./ (π * dosef)

    phi_val = repeat(phi_ph_val, 1, length(epsilon)) .+ transpose(phi_c_val)

    # Δ is derived from φ and Z in the solver after mixing; hand back φ here.
    return Zval, shift_val, phi_val
end


"""
    realEliashbergEq(mu_star, beta, znormip, phiip, shiftip, ws, w_static, w_static_chi, dosef, epsilon, dos, fermi_level)

Real axis Eliashberg equations in the vDOS+μ approximation. φ (= Δ·Z) is the primary
superconducting variable that is handed in and returned; Δ is derived from φ and Z in the solver.
"""
function realEliashbergEq(mu_star::Float64, beta::Float64, znormip::Vector{ComplexF64}, phiip::Vector{ComplexF64},
                          shiftip::Vector{ComplexF64}, ws::WPrimeWorkspace, w_static::AbstractVector,
                          w_static_chi, dosef::Float64, epsilon::Vector{Float64}, dos::Vector{Float64},
                          fermi_level::Float64, electronic_spec::Tuple)

    w_prime = ws.wp_full

    # interpolate Z,χ,ϕ onto ω'-integration grid (φ is handed in directly)
    Z_itp = linear_interpolation(w_static, znormip, extrapolation_bc=Flat())
    phi_itp = linear_interpolation(w_static, phiip, extrapolation_bc=Flat())
    shift_itp = linear_interpolation(w_static_chi, shiftip, extrapolation_bc=Flat())

    Z_ongrid = Z_itp.(w_prime)
    phi_ongrid = phi_itp.(w_prime)
    shift_ongrid = shift_itp.(w_prime) .- fermi_level

    # ------------- ε-integration ------------- #
    integrands = epsilon_helpers(electronic_spec, Z_ongrid, phi_ongrid, shift_ongrid, w_prime)


    # Z-integrand must be positive (causality), flip the others accordingly
    FLIP = ifelse.(integrands[1] .<= 0, 1, -1)
    g_z = -abs.(integrands[1])
    g_phi = -(FLIP .* integrands[2])
    g_chi = -(FLIP .* integrands[3])

    # ------------- ω'-integration ------------- #
    Iz = kernel_omega_integral(ws.Km_head, ws.Km_tail, ws, g_z)
    Iphi = kernel_omega_integral(ws.Kp_head, ws.Kp_tail, ws, g_phi)
    Ichi = kernel_omega_integral(ws.Kpchi_head, ws.Kpchi_tail, ws, g_chi)

    # μ* Coulomb term (scalar, same for all ω)
    coulomb = wprime_trapz(ws, g_phi .* tanh.(beta .* w_prime ./ 2))

    Zval = 1 .+ Iz ./ (w_static .* (π * dosef))
    phi_val = (Iphi .- mu_star .* coulomb) ./ (π * dosef)
    shift_val = -Ichi ./ (π * dosef)

    # Δ is derived from φ and Z in the solver after mixing; hand back φ here.
    return Zval, shift_val, phi_val
end



"""
    realEliashbergEq()

real axis Eliashberg equations in the vDOS+W approximation - OLD VERSION -
"""
function realEliashbergEq(beta::Float64, znormip::Vector{ComplexF64}, phiphip::Matrix{ComplexF64}, shiftip::Vector{ComplexF64},
                          Kp_func::AbstractInterpolation, Km_func::AbstractInterpolation, w_prime::Vector{Float64}, w_static::StepRangeLen,
                          w_static_chi, dosef::Float64, epsilon::Vector{Float64}, dos::Vector{Float64}, Weep::Matrix{Float64},
                          idx_ef::Int64, fermi_level::Float64, wgCoulomb::Number, electronic_spec::Tuple)

    if size(Weep, 1) != length(epsilon) || size(Weep, 2) != length(epsilon)
        error("The W(ε,ε′) matrix must be defined on the same energy grid as the DOS for real-axis vDOS+W calculations.")
    end

    # interpolate Z,χ,ϕ onto ω'-integration grid
    Z_itp = linear_interpolation(w_static, znormip, extrapolation_bc=Flat())
    shift_itp = linear_interpolation(w_static_chi, shiftip, extrapolation_bc=Flat())

    Z_ongrid = Z_itp.(w_prime)
    shift_ongrid = shift_itp.(w_prime) .- fermi_level
    phi_ongrid = interpolate_phi_matrix(w_static, phiphip, w_prime)

    # ------------- ε-integration ------------- #
    integrands, coulomb_spectral = eval_spectral_and_coulomb_vDOS_W(electronic_spec, epsilon, dos, Weep, Z_ongrid, phi_ongrid, shift_ongrid, w_prime)

    # ------------- Ω-integration ------------- #
    Kernel_minus, Kernel_plus = evaluate_Kernels(w_static, w_prime, Km_func, Kp_func)
    Kernel_plus_chi = evaluate_Kernels(w_static_chi, w_prime, Kp_func)
    z_integrand = -abs.(transpose(integrands[1])) .* Kernel_minus
    FLIP = ifelse.(-abs.(integrands[1]) .== integrands[1], 1, -1)
    phi_ph_integrand = -FLIP .* integrands[2] .* transpose(Kernel_plus)
    shift_integrand = -transpose(FLIP .* integrands[3]) .* Kernel_plus_chi

    # Coulomb term: integrate N(ε')W(ε,ε') with the same piecewise-linear
    # spectral quadrature, rather than multiplying by an interval-averaged W.
    coulomb_spectral = coulomb_spectral .* transpose(FLIP)
    coulomb_integrand = wgCoulomb .* dosef .* coulomb_spectral .* transpose(tanh.(beta .* w_prime ./ 2))

    # ------------- ω'-integration ------------- #
    Zval = 1 .+ 1 ./(w_static * π * dosef) .* trapz(w_prime, z_integrand)
    phi_ph_val = 1 ./(π * dosef) .* trapz(w_prime, transpose(phi_ph_integrand))
    phi_c_val = 1 ./(π * dosef) .* trapz(w_prime, coulomb_integrand)
    phi_val = repeat((phi_ph_val), 1, length(epsilon)) .+ transpose(phi_c_val)
    shift_val = -1 ./(π * dosef) .* trapz(w_prime, shift_integrand)

    delta_val = phi_val[:, idx_ef] ./ Zval


    return Zval, delta_val, shift_val, phi_val
end


"""
     realEliashbergEq(mu_star, beta, deltaip, ws, w_static)

Real axis Eliashberg equations in cDOS+μ approximation
"""
function realEliashbergEq(mu_star::Float64, beta::Float64, deltaip::Vector{ComplexF64},
                          ws::WPrimeWorkspace, w_static::AbstractVector)

    Delta_func = linear_interpolation(w_static, deltaip, extrapolation_bc=Flat())

    w_prime = ws.wp_full
    Delta_ongrid = Delta_func.(w_prime)
    sqrt_eval = @. sqrt(w_prime^2 - Delta_ongrid^2)

    g_z = @. real(w_prime / sqrt_eval)          # Z   integrand density (∝ quasiparticle DOS)
    g_phi = @. real(Delta_ongrid / sqrt_eval)   # Δ·Z integrand density

    # Ω & ω'-integration via the precomputed head/tail kernel blocks
    Iz = kernel_omega_integral(ws.Km_head, ws.Km_tail, ws, g_z)
    Iphi = kernel_omega_integral(ws.Kp_head, ws.Kp_tail, ws, g_phi)

    # μ* Coulomb term (scalar, same for all ω)
    coulomb = wprime_trapz(ws, g_phi .* tanh.(beta .* w_prime ./ 2))

    Zval = 1 .- Iz ./ w_static
    Delta = (Iphi .- mu_star .* coulomb) ./ Zval

    return Zval, Delta
end



###################################################
# ------------------- Helpers ------------------- #
###################################################
"""
    evaluate_Kernels(A, B, K_func)

Evaluation of the kernel at K(A,B)
"""
function evaluate_Kernels(A, B, K_func)
    N = length(A)
    M = length(B)
    K_vals = Matrix{ComplexF64}(undef, N, M)
    @inbounds for i in eachindex(B)
        b = B[i]   
        for j in eachindex(A)
            a = A[j]
            K_vals[j,i] = K_func(a, b)
        end
    end
    return K_vals
end

"""
    evaluate_Kernels(A, B, K_func)

Evaluation of the kernel at K(A,B)
"""
function evaluate_Kernels(A, B, K_func, K2_func)
    N = length(A)
    M = length(B)
    K_vals = Matrix{ComplexF64}(undef, N, M)
    K2_vals = Matrix{ComplexF64}(undef, N, M)
    @inbounds for i in eachindex(B)
        b = B[i]   
        for j in eachindex(A)
            a = A[j]
            K_vals[j,i] = K_func(a, b)
            K2_vals[j,i] = K2_func(a, b)
        end
    end
    return K_vals, K2_vals
end



function epsilon_helpers(electronic_spec, Z_ongrid, phi_ongrid, shift_ongrid, w_prime; idx_skip = nothing)
    eps_1, eps_2, dos_1, dos_2, dε, ddos = electronic_spec

    # assuming a linear form of the dos
    M1 = transpose((dos_2 - dos_1)./(eps_2 - eps_1))
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

    integrands = [Float64[] for _ in 1:3]
    integrands = Vector{Vector{Float64}}(undef, 3)
    for (idx, g) in enumerate((w_prime .* Z_ongrid, phi_ongrid, shift_ongrid))
        if ~isnothing(idx_skip) && any(idx .== idx_skip)
            continue
        end

        integrands[idx] = eval_spectral_integrals_fused(g, lorentzian_moments_p, lorentzian_moments_m, ε_p, Rplus, Rminus, Iplus, Iminus, idx == 3)
    end
                        

    return integrands
end

function interpolate_phi_matrix(w_static, phi::Matrix{ComplexF64}, w_prime::Vector{Float64})
    phi_ongrid = Matrix{ComplexF64}(undef, length(w_prime), size(phi, 2))
    @inbounds for iε in axes(phi, 2)
        phi_itp = linear_interpolation(w_static, phi[:, iε], extrapolation_bc=Flat())
        phi_ongrid[:, iε] = phi_itp.(w_prime)
    end
    return phi_ongrid
end


"""

Compute the Lorentzian integrals --> I should be able to sum here already over ε
"""
function Lorentzians_integrals_interval_grid(x0, x1, A::Vector{Float64}, B)

    N = length(x0)
    M = length(A)

    Int0 = Matrix{Float64}(undef, M, N)
    Int1 = Matrix{Float64}(undef, M, N)
    Int2 = Matrix{Float64}(undef, M, N)
    Int3 = Matrix{Float64}(undef, M, N)

    @inbounds for jε in 1:N
        xl = x0[jε]
        xr = x1[jε]

        @simd for iw in 1:M
            a = A[iw]
            b = B[iw]
            b2 = b^2
            a2 = a^2

            i0l = atan((a + xl) / b) / b
            i0r = atan((a + xr) / b) / b

            h_l = log((a + xl)^2 + b2)
            h_r = log((a + xr)^2 + b2)

            i1l = 0.5 * h_l - a * i0l
            i1r = 0.5 * h_r - a * i0r

            i2l = -a * h_l + (a2 - b2) * i0l + xl
            i2r = -a * h_r + (a2 - b2) * i0r + xr

            i3l = i0l * (3 * a * b2 - a * a2) + 0.5 * (3 * a2 - b2) * h_l + 0.5 * xl * (xl - 4 * a)
            i3r = i0r * (3 * a * b2 - a * a2) + 0.5 * (3 * a2 - b2) * h_r + 0.5 * xr * (xr - 4 * a)

            Int0[iw, jε] = i0r - i0l
            Int1[iw, jε] = i1r - i1l
            Int2[iw, jε] = i2r - i2l
            Int3[iw, jε] = i3r - i3l
        end
    end

    return Int0, Int1, Int2, Int3
end

function Lorentzians_integrals_interval_grid(x0, x1, A::Matrix{Float64}, B)

    Int0 = Matrix{Float64}(undef, size(A))
    Int1 = Matrix{Float64}(undef, size(A))
    Int2 = Matrix{Float64}(undef, size(A))
    Int3 = Matrix{Float64}(undef, size(A))

    @inbounds for jε in axes(A, 2)
        xl = x0[jε]
        xr = x1[jε]

        @simd for iw in axes(A, 1)
            a = A[iw, jε]
            b = B[iw, jε]
            b2 = b^2
            a2 = a^2

            i0l = atan((a + xl) / b) / b
            i0r = atan((a + xr) / b) / b

            h_l = log((a + xl)^2 + b2)
            h_r = log((a + xr)^2 + b2)

            i1l = 0.5 * h_l - a * i0l
            i1r = 0.5 * h_r - a * i0r

            i2l = -a * h_l + (a2 - b2) * i0l + xl
            i2r = -a * h_r + (a2 - b2) * i0r + xr

            i3l = i0l * (3 * a * b2 - a * a2) + 0.5 * (3 * a2 - b2) * h_l + 0.5 * xl * (xl - 4 * a)
            i3r = i0r * (3 * a * b2 - a * a2) + 0.5 * (3 * a2 - b2) * h_r + 0.5 * xr * (xr - 4 * a)

            Int0[iw, jε] = i0r - i0l
            Int1[iw, jε] = i1r - i1l
            Int2[iw, jε] = i2r - i2l
            Int3[iw, jε] = i3r - i3l
        end
    end

    return Int0, Int1, Int2, Int3
end

function epsilon_helpers_vDOSW(electronic_spec, epsilon, dos, wZ, phi_ongrid, shift, w_prime)
    eps_1, eps_2, dos_1, dos_2, dε, ddos = electronic_spec

    M1 = ddos ./ dε
    M0 = dos_1 .- eps_1 .* M1

    Φ1 = @views (phi_ongrid[:, 2:end] .- phi_ongrid[:,1:end-1]) ./ dε
    Φ0 = @views phi_ongrid[:, 1:end-1] .- eps_1 .* Φ1
    scale = 1 .+ Φ1.^2

    S = (shift .+ Φ0 .* Φ1) ./ scale
    P = sqrt.((wZ.^2 .- shift.^2 .- Φ0.^2) ./ scale .+ S.^2)

    Rplus = real.(S .+ P)
    Rminus = real.(S .- P)
    Iplus = imag.(S .+ P)
    Iminus = imag.(S .- P)

    I0_pp, I1_pp, I2_pp, I3_pp = Lorentzians_integrals_interval_grid(eps_1, eps_2, Rplus, Iplus)
    I0_mm, I1_mm, I2_mm, I3_mm = Lorentzians_integrals_interval_grid(eps_1, eps_2, Rminus, Iminus)
    
    return M0, M1, P, S, Φ0, Φ1, scale, Rplus, Rminus, Iplus, Iminus,
           I0_pp, I1_pp, I2_pp, I3_pp, I0_mm, I1_mm, I2_mm, I3_mm
end

function eval_spectral_and_coulomb_vDOS_W(electronic_spec, epsilon, dos, Weep, Z_ongrid, phi_ongrid, shift_ongrid, w_prime)
    eps_1, eps_2, dos_1, dos_2, dε, ddos = electronic_spec

    # necessary to reshape?
    wZ = reshape(w_prime .* Z_ongrid, :, 1)
    shift = reshape(shift_ongrid, :, 1)

    M0, M1, P, _, Φ0, Φ1, scale, Rplus, Rminus, Iplus, Iminus,
    I0_pp, I1_pp, I2_pp, I3_pp, I0_mm, I1_mm, I2_mm, I3_mm =
        epsilon_helpers_vDOSW(electronic_spec, epsilon, dos, wZ, phi_ongrid, shift, w_prime)

    inv_scale = 1 ./ scale
    
    z_integrand = eval_spectral_integrals(wZ .* inv_scale, M0, M1, P, Rplus, Rminus, Iplus, Iminus,
                                          I0_pp, I0_mm, I1_pp, I1_mm, I2_pp, I2_mm, I3_pp, I3_mm)
    phi_integrand = eval_spectral_integrals(Φ0 .* inv_scale, M0, M1, P, Rplus, Rminus, Iplus, Iminus,
                                            I0_pp, I0_mm, I1_pp, I1_mm, I2_pp, I2_mm, I3_pp, I3_mm, true;
                                            g_slope=Φ1 .* inv_scale)
    chi_integrand = eval_spectral_integrals(shift .* inv_scale, M0, M1, P, Rplus, Rminus, Iplus, Iminus,
                                            I0_pp, I0_mm, I1_pp, I1_mm, I2_pp, I2_mm, I3_pp, I3_mm, true;
                                            g_slope=inv_scale)

    NW0 = Weep[:, 1:end-1] .* dos_1
    NW1 = Weep[:, 2:end] .* dos_2
    M1W = (NW1 .- NW0) ./ dε
    M0W = NW0 .- eps_1 .* M1W

    C0, C1 = eval_spectral_integral_coefficients(Φ0 .* inv_scale, P, Rplus, Rminus, Iplus, Iminus,
                                                 I0_pp, I0_mm, I1_pp, I1_mm, I2_pp, I2_mm, I3_pp, I3_mm;
                                                 g_slope=Φ1 .* inv_scale)
    

    coulomb_spectral = Matrix{Float64}(undef, size(M0W, 1), size(C0, 1))
    mul!(coulomb_spectral, M0W, transpose(C0))
    mul!(coulomb_spectral, M1W, transpose(C1), 1.0, 1.0)

    return [z_integrand, phi_integrand, chi_integrand], coulomb_spectral
end


"""
    eval_speactral_integrals(g, I0_pp,  shift_int)
"""
function eval_spectral_integrals(g, M0, M1, ε_p, Rplus, Rminus, Iplus, Iminus, I0_pp, I0_mm, I1_pp, I1_mm, I2_pp, I2_mm, I3_pp, I3_mm, shift_int::Bool=false; g_slope=nothing)
    slope = if shift_int && isnothing(g_slope)
        one.(ε_p)
    else
        g_slope
    end
    C0, C1 = eval_spectral_integral_coefficients(g, ε_p, Rplus, Rminus, Iplus, Iminus,
                                                 I0_pp, I0_mm, I1_pp, I1_mm, I2_pp, I2_mm, I3_pp, I3_mm;
                                                 g_slope=slope)
    

    return vec(sum(M0 .* C0 .+ M1 .* C1, dims=2))
end


"""
    lorentzian_moments_weighted(x0, x1, A, B, M0, M1)

Evaluate lorentzian integrals together with weighting factors M0/M1
Allows to sum directly over epsilon
Returns M0_J^(k) in my notation
A, B are R+/I+ or R-/I-
"""
function lorentzian_moments_weighted(x0, x1, A::AbstractVector, B::AbstractVector, M0, M1)
    xl = vec(x0); xr = vec(x1); m0 = vec(M0); m1 = vec(M1)
    Nε = length(xl)
    M  = length(A)
    M0_I0 = zeros(M); M0_I1 = zeros(M); M0_I2 = zeros(M); M0_I3 = zeros(M)
    M1_I0 = zeros(M); M1_I1 = zeros(M); M1_I2 = zeros(M); M1_I3 = zeros(M)
    @inbounds for jε in 1:Nε
        l = xl[jε]; r = xr[jε]; w0 = m0[jε]; w1 = m1[jε]
        @simd for iw in 1:M
            a = A[iw]; b = B[iw]
            b2 = b^2; a2 = a^2
            i0l = atan((a + l) / b) / b;  i0r = atan((a + r) / b) / b
            hl  = log((a + l)^2 + b2);    hr  = log((a + r)^2 + b2)
            i1l = 0.5 * hl - a * i0l;     i1r = 0.5 * hr - a * i0r
            i2l = -a * hl + (a2 - b2) * i0l + l;  i2r = -a * hr + (a2 - b2) * i0r + r
            i3l = i0l * (3 * a * b2 - a * a2) + 0.5 * (3 * a2 - b2) * hl + 0.5 * l * (l - 4 * a)
            i3r = i0r * (3 * a * b2 - a * a2) + 0.5 * (3 * a2 - b2) * hr + 0.5 * r * (r - 4 * a)
            d0 = i0r - i0l; d1 = i1r - i1l; d2 = i2r - i2l; d3 = i3r - i3l
            M0_I0[iw] += w0 * d0; M0_I1[iw] += w0 * d1; M0_I2[iw] += w0 * d2; M0_I3[iw] += w0 * d3
            M1_I0[iw] += w1 * d0; M1_I1[iw] += w1 * d1; M1_I2[iw] += w1 * d2; M1_I3[iw] += w1 * d3
        end
    end

    return (M0_I0, M0_I1, M0_I2, M0_I3, M1_I0, M1_I1, M1_I2, M1_I3)
end


"""
    eval_spectral_integrals_fused(g, M0, M1, ε_p, Rplus, Rminus, Iplus, Iminus, eps_1, eps_2, shift_int; g_slope)

Evaluate spectral integrals for wZ, delta, chi based on the precomputed lorentzian integrals
"""
function eval_spectral_integrals_fused(g, lorentzian_moments_p, lorentzian_moments_m, ε_p, Rplus, Rminus, Iplus, Iminus, shift_int::Bool=false; g_slope=nothing)
    slope = (shift_int && isnothing(g_slope)) ? one.(ε_p) : g_slope

    M0_I0_pp, M0_I1_pp, M0_I2_pp, M0_I3_pp, M1_I0_pp, M1_I1_pp, M1_I2_pp, M1_I3_pp = lorentzian_moments_p
    M0_I0_mm, M0_I1_mm, M0_I2_mm, M0_I3_mm, M1_I0_mm, M1_I1_mm, M1_I2_mm, M1_I3_mm = lorentzian_moments_m

    # spectral coefficients Rg, Ig, R+-,I+-
    C0s, _ = eval_spectral_integral_coefficients(g, ε_p, Rplus, Rminus, Iplus, Iminus,
                                                 M0_I0_pp, M0_I0_mm, M0_I1_pp, M0_I1_mm, M0_I2_pp, M0_I2_mm, M0_I3_pp, M0_I3_mm; g_slope=slope)
    _, C1s = eval_spectral_integral_coefficients(g, ε_p, Rplus, Rminus, Iplus, Iminus,
                                                 M1_I0_pp, M1_I0_mm, M1_I1_pp, M1_I1_mm, M1_I2_pp, M1_I2_mm, M1_I3_pp, M1_I3_mm; g_slope=slope)
    return C0s .+ C1s
end


function eval_spectral_integral_coefficients(g, ε_p, Rplus, Rminus, Iplus, Iminus, I0_pp, I0_mm, I1_pp, I1_mm, I2_pp, I2_mm, I3_pp, I3_mm; g_slope=nothing)
    Rg = real(g./ (2*ε_p))
    Ig = imag(g./ (2*ε_p))

    C0 = @. Ig * (I1_pp - I1_mm)
    C1 = @. Ig * (I2_pp - I2_mm)
    @. C0 += Ig * (Rplus * I0_pp - Rminus * I0_mm)
    @. C1 += Ig * (Rplus * I1_pp - Rminus * I1_mm)
    @. C0 += Rg * (-Iplus * I0_pp + Iminus * I0_mm)
    @. C1 += Rg * (-Iplus * I1_pp + Iminus * I1_mm)

    if !isnothing(g_slope)
        # Add the ε-dependent part of an affine numerator g(ε)=g+g_slope*ε.
        # The vDOS+μ χ integral corresponds to g_slope=1.
        Rg = real(g_slope ./ (2 * ε_p))
        Ig = imag(g_slope ./ (2 * ε_p))

        @. C0 += Ig * (I2_pp - I2_mm)
        @. C1 += Ig * (I3_pp - I3_mm)
        @. C0 += Ig * (Rplus * I1_pp - Rminus * I1_mm)
        @. C1 += Ig * (Rplus * I2_pp - Rminus * I2_mm)
        @. C0 += -Rg * (Iplus * I1_pp - Iminus * I1_mm)
        @. C1 += -Rg * (Iplus * I2_pp - Iminus * I2_mm)
    end

    return C0, C1
end
