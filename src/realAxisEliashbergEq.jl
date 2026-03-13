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

real axis Eliashberg equations in vDOS+μ approximation
"""
function realEliashbergEq(mu_star::Float64, beta::Float64, znormip::Vector{ComplexF64}, deltaip::Vector{ComplexF64}, shiftip::Vector{ComplexF64},
                          Kp_func::AbstractInterpolation, Km_func::AbstractInterpolation, w_dynam::Vector{Float64}, w_static::StepRangeLen, 
                          w_cut::Float64, dosef::Float64, epsilon::Vector{Float64}, dos::Vector{Float64}, fermi_level::Float64)

    # delta/Z
    phiphip = deltaip ./ znormip
    
    # shift dos grid to fermi level
    epsilon = epsilon .- fermi_level   
    
    # Find root of Θ(ω')
    theta = @. w_static'^2*znormip'^2 - (epsilon+shiftip')^2 - phiphip'^2
    #theta_func = linear_interpolation(w_static, theta, extrapolation_bc=Flat())
    itp = interpolate(theta, (NoInterp(), BSpline(Linear())))
    theta_func = extrapolate(scale(itp, 1:length(epsilon), w_static), Flat())
    gap0 = real(phiphip[1]/znormip[1])

    root_eq = x -> real(theta_func(:,x))    # !!!!!!! Does not work this way. Try to perform the epsilon integration first and then interpolate !!!!!!!!!!!1
    root = find_zero(root_eq, gap0)

    # shift chebyshev grid to pole of Θ(ω') 
    w_prime = w_dynam .+ root   
    w_pm = ifelse.((w_prime .< w_cut) .& (w_prime .> 0), w_prime, 0.0)
    
    
    # ---------- ε-integration ---------- # 
    # using a linear interpolation of the dos
    M1 = (dos[2:end] - dos[1:end-1])/(epsilon[2:end] - epsilon[1:end-1])
    M0 = dos[1:end-1] - epsilon[1:end-1] * M1
    
    # Helper functions for the ε-integration
    ε_p = sqrt(w_pm.^2*znormip.^2 - phiphip.^2)
    Rplus = real(shiftip .+ ε_p)
    Rminus = real(shiftip .- ε_p)
    Iplus = imag(shiftip .+ ε_p)
    Iminus = imag(shiftip .- ε_p)

    # Lorentzian integrals
    # pp --> ++ signature (Rplus, Iplus)
    I0_pp, I1_pp, I2_pp = Lorentzians_integrals(epsilon, Rplus, Iplus)
    I0_mm, I1_mm, I2_mm = Lorentzians_integrals(epsilon, Rminus, Iminus)
    
    # ω' integrands
    integrands = Vector{Vector{ComplexF64}}(undef, 3)
    for (idx, g) in enumerate((znormip, phiphip, shiftip))
        Rg = real(g/2/ε_p)
        Ig = imag(g/2/ε_p)

        integrands[idx] += @. M0 * Ig *(I1_pp - I1_mm) + M1*Ig * (I2_pp -I2_mm)
        integrands[idx] += @. M0 *Ig *(Rplus*Io_pp - Rminus*I0_mm) + M1*Ig *(Rplus*I1_pp - Rminus*I1_mm)
        
        # shift contains a ε-dependence in the numerator
        if idx == 3
            integrands[idx] += @. - M0 * Rg *(Iplus*I1_pp - Iminus*I1_mm) - M1*Rg *(Iplus*I2_pp - Iminus*I2_mm)
        else
            integrands[idx] += @. - M0*Rg*(Iplus*I0_pp -Iminus*I0_mm) - M1*Rg*(Iplus*I1_pp - Iminus*I1_mm)
        end
    end

    # ---------- Ω-integration ---------- #
    # evaluate K(ω,ω')
    # !!!!!!!!!!!!!!!!! Chekc if Km/Kp is correct in Z,chi, Delta !!!!!!!!!!!!!!!!!
    evalKernel = evaluate_Kernels(w_static, w_pm, Km_func)
    z_integrand = integrands[1] .* evalKernel'

    evalKernel = evaluate_Kernels(w_static, w_pm, Kp_func)
    phi_integrand = integrands[2] .* (evalKernel' + 1/2 *mu_star * tanh.(beta * w_pm / 2))
    shift_integrand = integrands[3] .* evalKernel'

    # ω'-integration
    Zval = 1 .+ 1/(w_static* π *dosef)*trapz(w_prime, z_integrand)
    phi_val = 1/(π *dosef)*trapz(w_prime, phi_integrand)
    shift_val = 1/(π *dosef)*trapz(w_prime, shift_integrand)

    delta_val = phi_val ./ Zval

    return Zval, delta_val, shift_val
    
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
    w_pm = ifelse.((w_prime .< w_cut) .& (w_prime .> 0), w_prime, 0.0)

    # delta on chebyshev grid
    Delta_func_eval = Delta_func.(w_pm)
    sqrt_eval = @. sqrt((w_pm^2 - Delta_func_eval^2))

    # Evaluate Z
    evalKernel = evaluate_Kernels(w_static, w_pm, Km_func)
    realPart = @. real(w_pm / sqrt_eval)
    integrand = @. realPart * evalKernel

    Zval = 1 .- trapz(w_prime, transpose(integrand)) ./ w_static

    # Evaluate Delta
    eval_real = @. real(Delta_func_eval / sqrt_eval)

    evalKernel = evaluate_Kernels(w_static, w_pm, Kp_func)
    int1 = @. eval_real * evalKernel
    int2 = @. mu_star *eval_real * tanh(beta * w_pm / 2)

    Delta = (trapz(w_prime, transpose(int1)) .- trapz(w_prime, transpose(int2))) ./ Zval

    return Zval, Delta

end



###################################################
# ------------------- Helpers ------------------- #   
###################################################
"""
    evaluate_Kernels(A, B, K_func)

Evaluation of the kernel at K(A,B)
"""
function evaluate_Kernels(A::StepRangeLen, B::Vector{Float64}, K_func)
    N = length(A)
    M = length(B)
    K_vals = Matrix{ComplexF64}(undef, N, M)
    @inbounds for i in eachindex(A) 
        a = A[i]
        for j in eachindex(B)
            b = B[j]
            K_vals[j,i] = K_func(a, b)
        end
    end
    return K_vals
end


""" 
    Lorentzians_integrals()
        
Analytical epxressions for the integrals over Lorentzians times a Polynomial, used
to evalaute the spectral integrals.
"""
function Lorentzians_integrals(x, A, B)

    # I1: 1/((x+A)^2 + B^2)
    I0 = arctan((A+x)/B)/B

    # I2: x/((x+A)^2 + B^2)
    I1h = log((A+x)^2 + B^2)
    I1 = 0.5*I1h - A*I0

    # I3: x^2/((x+A)^2 + B^2)
    I2 = -A*I1h + (A^2-B^2)*I0 + x

    return I0, I1, I2

end