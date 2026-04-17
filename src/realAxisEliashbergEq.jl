"""
    file containing the real axis eliashberg equations
"""


"""

real axis Eliashberg equations in the vDOS+W approximation
"""
function realEliashbergEq(itemp::Number, a2f_omega::StepRangeLen{Float64}, a2f::Vector{Float64}, dosef::Float64, ndos::Int64, dos_en::Vector{Float64}, dos::Vector{Float64}, 
                            Weep::Matrix{Float64}, znormip::Vector{ComplexF64}, phiphip::Vector{ComplexF64}, phicip::Vector{ComplexF64}, shiftip::Vector{ComplexF64}, wgCoulomb::Float64)
    """
    vDOS + W equations:
    - itemp: current temperature
    _ a2f_omega: alpha^2F energy grid
    - a2f: alpha^2F values
    - dosef: dos at fermi energy
    - ndos: length of dos vector
    - dos_en: energy grid dos
    - dos: dos values
    - Weep: W(ε,ε')-values, matrix
    - znormip: Z(ω) previous iteration
    - phiphip: ϕ_ph(ω) previous iteration, phonon part
    - phicip: ϕ_c(ω) previous iteration, coulomb part
    - shiftip: χ(ω) previous iterations
    - wgCoulomb: weigth coulomb interaction, damping during first iterations
    """


    # dummy variables
    znormi = ones(ComplexF64, length(znormip))
    phiphi = ones(ComplexF64, length(phiphip))
    phici = ones(ComplexF64, length(phicip))
    shifti = ones(ComplexF64, length(shiftip))


    data = Vector{Vector{ComplexF64}}([znormi, phiphi, phici, shifti])
    return data 

end


"""
    realEliashbergEq()

real axis Eliashberg equations in vDOS+μ approximation
"""
function realEliashbergEq(mu_star::Float64, beta::Float64, znormip::Vector{ComplexF64}, deltaip::Vector{ComplexF64}, shiftip::Vector{ComplexF64},
                          Kp_func::AbstractInterpolation, Km_func::AbstractInterpolation, w_prime::Vector{Float64}, w_static::StepRangeLen, 
                          w_static_chi, dosef::Float64, epsilon::Vector{Float64}, dos::Vector{Float64}, fermi_level::Float64, i_it)

    # delta/Z
    phiphip = deltaip .* znormip
    
    # shift dos grid to fermi level
    #epsilon = epsilon .- fermi_level 

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
        
        integrands[idx] = eval_spectral_integrals(g, M0, M1, ε_p, Rplus, Rminus, Iplus, Iminus, I0_pp, I0_mm, I1_pp, I1_mm, I2_pp, I2_mm, I3_pp, I3_mm, shift_int, python_style=true, epsilon=epsilon)
        #println(maximum(integrands[idx]))
    end

    #println((integrands[1][1:5]))
 
    # ------------- Ω-integration ------------- #
    # evaluate K(ω,ω')
    # !!!!!!!!!!!!!!!!! Check if Km/Kp is correct in Z,chi, Delta !!!!!!!!!!!!!!!!!
    Kernel_minus = evaluate_Kernels(w_static, w_prime, Km_func)
    Kernel_plus = evaluate_Kernels(w_static, w_prime, Kp_func)
    Kernel_plus_chi = evaluate_Kernels(w_static_chi, w_prime, Kp_func)


    z_integrand = -abs.(transpose(integrands[1])) .* Kernel_minus    # Z-integrand must be positive (causality)
    FLIP = ifelse.(-abs.(integrands[1]) .== integrands[1], 1, -1)    # enforce causality through sign flip
    phi_integrand = -FLIP.*integrands[2] .* (transpose(Kernel_plus) .+ mu_star * tanh.(beta .* w_prime ./ 2))
    shift_integrand = -transpose(FLIP.*integrands[3]) .* Kernel_plus_chi

    # println("sum ", (sum(trapz(w_prime, Kernel_minus))))

    # if i_it <3
    #     plot(w_prime, (-abs.((integrands[1]))), label="real")
    #     savefig("Z_integrand_$i_it.png")

    #     plot(w_prime, real(Kernel_minus[11,:]), label="real")
    #     plot!(w_prime, imag(Kernel_minus[11,:]), label="imag")
    #     savefig("Kminus_$i_it.png")

    #     plot(w_prime, real(Kernel_minus[320,:]), label="real")
    #     plot!(w_prime, imag(Kernel_minus[320,:]), label="imag")
    #     savefig("Kminus_320_$i_it.png")

    #     plot(w_prime, ((real(z_integrand[1,:]))), label="real")
    #     plot!(w_prime, ((imag(z_integrand[1,:]))), label="imag")
    #     savefig("ZK_integrand_$i_it.png")

    #     plot(w_prime, ((integrands[2])), label="real")
    #     savefig("phi_integrand_$i_it.png")

    #     plot(w_prime, ((integrands[3])), label="real")
    #     savefig("shift_integrand_$i_it.png")

    #     plot(w_static, real(trapz(w_prime, z_integrand)), label="real")
    #     plot!(w_static, imag(trapz(w_prime, z_integrand)), label="imag")
    #     savefig("Zintegrated_$i_it.png")

    #     plot(w_static, real(1 ./(w_static* π *dosef).*trapz(w_prime, z_integrand)), label="real")
    #     plot!(w_static, imag(1 ./(w_static* π *dosef) .*trapz(w_prime, z_integrand)), label="imag")
    #     savefig("AllZintegrated_$i_it.png")
    # end


    # ------------- ω'-integration ------------- #
    # pole is automatically at 0 after ε-integration
    Zval = 1 .+ 1 ./(w_static* π *dosef) .*trapz(w_prime, z_integrand)  
    phi_val = 1 ./(π *dosef) .*trapz(w_prime, transpose(phi_integrand))
    shift_val = -1 ./(π *dosef) .*trapz(w_prime, shift_integrand)
    #shift_val .= 0

    delta_val = phi_val ./ Zval

    return Zval, delta_val, shift_val, FLIP
    
end



"""
     realEliashbergEq(mu_star, beta, Delta_func_eval, sqrt_eval, Kp_func, Km_func, w_prime, w_max, W_pm, w_static, W_sm, W_prime) 

Real axis Eliashberg equations in cDOS+μ approximation
"""
function realEliashbergEq(mu_star::Float64, beta::Float64, deltaip::Vector{ComplexF64}, Kp_func::AbstractInterpolation, 
                          Km_func::AbstractInterpolation, w_dynam::Vector{Float64}, w_static::StepRangeLen, w_cut::Float64) 

    # w_static: grid of the static part of the kernels, on which Z and Δ are evaluated
    # w_dynam:  chebyshev grid shifted to the pole of Θ(ω') to capture the singularity,
    #           used for the integration

    Delta_func = linear_interpolation(w_static, deltaip, extrapolation_bc=Flat())
    gap0 = real(deltaip[1])

    root_eq = x -> x - real(Delta_func(x))
    root = find_zero(root_eq, gap0)

    # update 
    w_prime = w_dynam .+ root    # shift chebyshev grid to pole of Θ(ω')
    # w_pm = ifelse.((w_prime .< w_cut) .& (w_prime .> 0), w_prime, 0.0)    # MIT version
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

    evalKernel = transpose(evaluate_Kernels( w_static, w_pm, Kp_func))
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
    @inbounds for i in eachindex(A)   
        a = A[i]
        for j in eachindex(B)
            b = B[j]
            K_vals[i,j] = K_func(a, b)
        end
    end
    return K_vals
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


"""
    eval_speactral_integrals(g, I0_pp,  shift_int)
"""
function eval_spectral_integrals(g, M0, M1, ε_p, Rplus, Rminus, Iplus, Iminus, I0_pp, I0_mm, I1_pp, I1_mm, I2_pp, I2_mm, I3_pp, I3_mm, shift_int::Bool=false; python_style::Bool=false, epsilon=nothing)
   
    Rg = real(g./ (2*ε_p))
    Ig = imag(g./ (2*ε_p))

    if python_style
        isnothing(epsilon) && error("epsilon must be supplied when python_style=true")

        E1 = transpose(epsilon[1:end-1])
        E2 = transpose(epsilon[2:end])

        function dos_int_terms(R, I, Rg, Ig)
            B = @. R*Ig - I*Rg
            log_diff = @. log((E2 + R)^2 + I^2) - log((E1 + R)^2 + I^2)
            atan_diff = @. atan((E2 + R)/I) - atan((E1 + R)/I)

            C_log = @. 0.5*(M0*Ig + M1*B) - M1*Ig*R
            C_atan = @. M0*B/I + M1*Ig*(R^2 - I^2)/I - R/I*(M0*Ig + M1*B)
            C_x = @. M1*Ig

            return @. C_log*log_diff + C_atan*atan_diff + C_x*(E2 - E1)
        end

        integrand_temp = dos_int_terms(Rplus, Iplus, Rg, Ig) .- dos_int_terms(Rminus, Iminus, Rg, Ig)

        if shift_int
            Rg_eps = real(1 ./ (2*ε_p))
            Ig_eps = imag(1 ./ (2*ε_p))

            function dos_int_eps_terms(R, I, Rg, Ig)
                B = @. R*Ig - I*Rg
                log_diff = @. log((E2 + R)^2 + I^2) - log((E1 + R)^2 + I^2)
                atan_diff = @. atan((E2 + R)/I) - atan((E1 + R)/I)

                C_log = @. -R*(M0*Ig + M1*B) + 0.5*M0*B + 0.5*M1*Ig*(3*R^2 - I^2)
                C_atan = @. M1*Rg*(I^2 - R^2) + 2*M1*Ig*R*I + M0*(Rg*R - Ig*I)
                C_x2 = @. 0.5*M1*Ig
                C_x = @. M0*Ig + M1*B

                return @. C_log*log_diff + C_atan*atan_diff + C_x2*(E2*(E2 - 4*R) - E1*(E1 - 4*R)) + C_x*(E2 - E1)
            end

            integrand_temp += dos_int_eps_terms(Rplus, Iplus, Rg_eps, Ig_eps) .- dos_int_eps_terms(Rminus, Iminus, Rg_eps, Ig_eps)
        end

    else

        integrand_temp =  @. M0 * Ig *(I1_pp - I1_mm) + M1*Ig * (I2_pp -I2_mm)
        integrand_temp += @. M0 *Ig *(Rplus*I0_pp - Rminus*I0_mm) + M1*Ig *(Rplus*I1_pp - Rminus*I1_mm)
        integrand_temp += @. M0*Rg*(-Iplus*I0_pp +Iminus*I0_mm) + M1*Rg*(-Iplus*I1_pp + Iminus*I1_mm)

        if shift_int
            # shift contains a ε-dependence in the numerator (Rg)
            Rg = real(1 ./ (2*ε_p))
            Ig = imag(1 ./ (2*ε_p))

            integrand_temp +=  @. M0 * Ig *(I2_pp - I2_mm) + M1*Ig * (I3_pp -I3_mm)
            integrand_temp += @. M0 *Ig *(Rplus*I1_pp - Rminus*I1_mm) + M1*Ig *(Rplus*I2_pp - Rminus*I2_mm)
            integrand_temp += @.  - M0 .*Rg  .*(Iplus.*I1_pp .- Iminus.*I1_mm) .- M1 .*Rg .*(Iplus.*I2_pp .- Iminus.*I2_mm)
        end

        
    end

    #println("AD: ", integrand_temp[1:5,1:5])

    return vec(sum(integrand_temp, dims=2))

end