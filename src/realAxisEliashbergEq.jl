"""
    file containing the real axis eliashberg equations
"""


"""

real axis Eliashberg equations in the vDOS+W approximation
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
    realEliashbergEq()

real axis Eliashberg equations in vDOS+μ approximation
"""
function realEliashbergEq(mu_star::Float64, beta::Float64, znormip::Vector{ComplexF64}, deltaip::Vector{ComplexF64}, shiftip::Vector{ComplexF64},
                          Kp_func::AbstractInterpolation, Km_func::AbstractInterpolation, w_prime::Vector{Float64}, w_static::StepRangeLen, 
                          w_static_chi, dosef::Float64, epsilon::Vector{Float64}, dos::Vector{Float64}, fermi_level::Float64)

    # delta/Z
    phiphip = deltaip .* znormip

    # interpolate Z,χ,ϕ onto ω'-integration grid
    Z_itp = linear_interpolation(w_static, znormip, extrapolation_bc=Flat())
    phi_itp = linear_interpolation(w_static, phiphip, extrapolation_bc=Flat())
    shift_itp = linear_interpolation(w_static_chi, shiftip, extrapolation_bc=Flat())

    Z_ongrid = Z_itp.(w_prime)
    phi_ongrid = phi_itp.(w_prime)
    shift_ongrid = shift_itp.(w_prime) .- fermi_level

    
    # ------------- ε-integration ------------- # 
    M0, M1, ε_p, Rplus, Rminus, Iplus, Iminus, I0_pp, I1_pp, I2_pp, I3_pp, I0_mm, I1_mm, I2_mm, I3_mm = epsilon_helpers(epsilon, dos, Z_ongrid, phi_ongrid, shift_ongrid, w_prime)
    
    # ω' integrands
    integrands = Vector{Vector{Float64}}(undef, 3)
    integrands = [zeros(size(w_prime)) for _ in 1:3]

    for (idx, g) in enumerate((w_prime.*Z_ongrid, phi_ongrid, shift_ongrid))

        shift_int = false
        if idx == 3
            shift_int = true    # extra terms in ε'-integration
        end
        
        integrands[idx] = eval_spectral_integrals(g, M0, M1, ε_p, Rplus, Rminus, Iplus, Iminus, I0_pp, I0_mm, I1_pp, I1_mm, I2_pp, I2_mm, I3_pp, I3_mm, shift_int)
    end

    # ------- testing ------- #
    pole_eqs = [
        x -> x.^2 .*imag(Z_itp(x)) .*real(Z_itp(x)) - imag(phi_itp(x)) .*real(phi_itp(x)),
        x -> x.^2 .*(real(Z_itp(x)).^2 - imag(Z_itp(x).^2)) - real(phi_itp(x)).^2 + imag(phi_itp(x)).^2
        ]
    poles = []
    for pole_eq in pole_eqs
        push!(poles, find_zero(pole_eq, gap0))
    end
    plot(w_prime, pole_eqs[1](w_primeas))
    savefig("pole_eq.png")
    error("A")
 
    # ------------- Ω-integration ------------- #
    # evaluate K(ω,ω')
    Kernel_minus = evaluate_Kernels(w_static, w_prime, Km_func)
    Kernel_plus = evaluate_Kernels(w_static, w_prime, Kp_func)
    Kernel_plus_chi = evaluate_Kernels(w_static_chi, w_prime, Kp_func)

    z_integrand = -abs.(transpose(integrands[1])) .* Kernel_minus    # Z-integrand must be positive (causality)
    FLIP = ifelse.(-abs.(integrands[1]) .== integrands[1], 1, -1)    # enforce causality through sign flip
    phi_integrand = -FLIP.*integrands[2] .* (transpose(Kernel_plus) .- mu_star * tanh.(beta .* w_prime ./ 2))
    shift_integrand = -transpose(FLIP.*integrands[3]) .* Kernel_plus_chi

    # ------------- ω'-integration ------------- #
    # pole is at 0 after ε-integration
    Zval = 1 .+ 1 ./(w_static* π *dosef) .*trapz(w_prime, z_integrand)  
    phi_val = 1 ./(π *dosef) .*trapz(w_prime, transpose(phi_integrand))
    shift_val = -1 ./(π *dosef) .*trapz(w_prime, shift_integrand)

    delta_val = phi_val ./ Zval

    return Zval, delta_val, shift_val
    
end



"""
     realEliashbergEq(mu_star, beta, Delta_func_eval, sqrt_eval, Kp_func, Km_func, w_prime, w_max, W_pm, w_static, W_sm, W_prime) 

Real axis Eliashberg equations in cDOS+μ approximation
"""
function realEliashbergEq(mu_star::Float64, beta::Float64, deltaip::Vector{ComplexF64}, Kp_func::AbstractInterpolation, 
                          Km_func::AbstractInterpolation, w_dynam::Vector{Float64}, w_static::StepRangeLen, w_cut::Float64, i_it) 

    # w_static: grid of the static part of the kernels, on which Z and Δ are evaluated
    # w_dynam:  chebyshev grid shifted to the pole of Θ(ω') to capture the singularity,
    #           used for the integration

    Delta_func = linear_interpolation(w_static, deltaip, extrapolation_bc=Flat())
    gap0 = real(deltaip[1])

    root_eq = x -> x - real(Delta_func(x))

    root = gap0
    try
        root = find_zero(root_eq, gap0)
    catch
        root = find_zero(root_eq, 20)
    end

    # update 
    w_prime = w_dynam .+ root    # shift chebyshev grid to pole of Θ(ω')
    w_pm = w_prime[(w_prime .< w_cut) .& (w_prime .> 0)]

    # delta on chebyshev grid
    Delta_func_eval = Delta_func.(w_pm)
    sqrt_eval = @. sqrt((w_pm^2 - Delta_func_eval^2))

    # Evaluate Z
    evalKernel = transpose(evaluate_Kernels(w_static, w_pm,  Km_func))
    realPart = @. real(w_pm / sqrt_eval)
    integrand = @. realPart * evalKernel

    Zval = 1 .- trapz(w_pm, transpose(integrand)) ./ w_static

    # Evaluate Delta
    eval_real = @. real(Delta_func_eval / sqrt_eval)

    evalKernel = transpose(evaluate_Kernels(w_static, w_pm, Kp_func))
    int1 = @. eval_real * evalKernel
    int2 = @. mu_star *eval_real * tanh(beta * w_pm / 2)

    Delta = (trapz(w_pm, transpose(int1)) .- trapz(w_pm, transpose(int2))) ./ Zval

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



""" 
    Lorentzians_integrals()
        
Analytical epxressions for the integrals over Lorentzians times a Polynomial, used
to evaluate the spectral integrals.
"""
function Lorentzians_integrals(x::Vector{Float64}, A::Vector{Float64}, B::Vector{Float64})

    x = transpose(x[:])

    # I0: 1/((x+A)^2 + B^2)
    I0 = @. atan((A+x)/B)/B

    # I1: x/((x+A)^2 + B^2)
    I1h = @. log((A+x)^2 + B^2)
    I1 = @. 0.5*I1h - A*I0

    # I2: x^2/((x+A)^2 + B^2)
    I2 = @. -A*I1h + (A^2-B^2)*I0 + x

    # I3: x^3/((x+A)^2 + B^2)
    I3 = @. I0*(3*A*B^2-A^3) + 0.5 *(3*A^2-B^2)*I1h + 0.5*x*(x-4*A)

    # Integration boundaries [ej, ej+1]
    Int0 = I0[:, 2:end] .- I0[:, 1:end-1]
    Int1 = I1[:, 2:end] .- I1[:, 1:end-1]
    Int2 = I2[:, 2:end] .- I2[:, 1:end-1]
    Int3 = I3[:, 2:end] .- I3[:, 1:end-1]

    return Int0, Int1, Int2, Int3
end


function epsilon_helpers(epsilon, dos, Z_ongrid, phi_ongrid, shift_ongrid, w_prime)
    # assuming a linear form of the dos
    M1 = transpose((dos[2:end] .- dos[1:end-1])./(epsilon[2:end] .- epsilon[1:end-1]))
    M0 = transpose(dos[1:end-1]) .- transpose(epsilon[1:end-1]) .* M1
    
    # Helper functions for the ε-integration
    ε_p = sqrt.(w_prime.^2 .*Z_ongrid.^2 .- phi_ongrid.^2)       
    Rplus = real.(shift_ongrid .+ ε_p)
    Rminus = real.(shift_ongrid .- ε_p)
    Iplus = imag.(shift_ongrid .+ ε_p)
    Iminus = imag.(shift_ongrid .- ε_p)

    # Lorentzian integrals
    # pp --> ++ signature (Rplus, Iplus)
    I0_pp, I1_pp, I2_pp, I3_pp = Lorentzians_integrals(epsilon, Rplus, Iplus)
    I0_mm, I1_mm, I2_mm, I3_mm = Lorentzians_integrals(epsilon, Rminus, Iminus)

    return M0, M1, ε_p, Rplus, Rminus, Iplus, Iminus, I0_pp, I1_pp, I2_pp, I3_pp, I0_mm, I1_mm, I2_mm, I3_mm
end

function interpolate_phi_matrix(w_static, phi::Matrix{ComplexF64}, w_prime::Vector{Float64})
    phi_ongrid = Matrix{ComplexF64}(undef, length(w_prime), size(phi, 2))
    @inbounds for iε in axes(phi, 2)
        phi_itp = linear_interpolation(w_static, phi[:, iε], extrapolation_bc=Flat())
        phi_ongrid[:, iε] = phi_itp.(w_prime)
    end
    return phi_ongrid
end

function Lorentzians_integrals_interval_grid(x0, x1, A, B)

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

# Slightly faster (maybe)
    # println("new")
    # @time begin
    #     nω = size(phi_ongrid, 1)
    #     nε = size(phi_ongrid, 2) - 1

    #     Φ1 = Matrix{ComplexF64}(undef, nω, nε)
    #     Φ0 = Matrix{ComplexF64}(undef, nω, nε)
    #     scale = Matrix{ComplexF64}(undef, nω, nε)
    #     S = Matrix{ComplexF64}(undef, nω, nε)
    #     P = Matrix{ComplexF64}(undef, nω, nε)

    #     Rplus = Matrix{Float64}(undef, nω, nε)
    #     Rminus = Matrix{Float64}(undef, nω, nε)
    #     Iplus = Matrix{Float64}(undef, nω, nε)
    #     Iminus = Matrix{Float64}(undef, nω, nε)

    #     @inbounds for jε in 1:nε
    #         e1 = eps_1[jε]
    #         de = dε[jε]

    #         @simd for iw in 1:nω
    #             phi1 = (phi_ongrid[iw, jε + 1] - phi_ongrid[iw, jε]) / de
    #             phi0 = phi_ongrid[iw, jε] - e1 * phi1
    #             sc = 1 + phi1^2

    #             s = (shift[iw] + phi0 * phi1) / sc
    #             p = sqrt((wZ[iw]^2 - shift[iw]^2 - phi0^2) / sc + s^2)

    #             Φ1[iw, jε] = phi1
    #             Φ0[iw, jε] = phi0
    #             scale[iw, jε] = sc
    #             S[iw, jε] = s
    #             P[iw, jε] = p

    #             sp = s + p
    #             sm = s - p
    #             Rplus[iw, jε] = real(sp)
    #             Rminus[iw, jε] = real(sm)
    #             Iplus[iw, jε] = imag(sp)
    #             Iminus[iw, jε] = imag(sm)
    #         end
    #     end
    # end

    # println(all((Rplus2 -Rplus) .< 1e-12))
    # println(all((Rminus2 - Rminus) .< 1e-12))
    # println(all((Iplus2 - Iplus) .< 1e-12))
    # println(all((Iminus2 - Iminus) .< 1e-12))


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
