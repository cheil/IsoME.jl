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
function realEliashbergEq(mu::Float64, beta::Float64, Delta_func_eval::Vector{ComplexF64}, sqrt_eval::Vector{ComplexF64}, 
                          Kp_func::AbstractInterpolation, Km_func::AbstractInterpolation, w_prime::Vector{Float64}, 
                          w_static::StepRangeLen, w_pm::Vector{Float64}, dosef::Float64, ndos::Int64, dos_en::Vector{Float64},
                          dos::Vector{Float64})

    # Evaluate Z
    evalKernel = evaluate_Kernels(w_static, w_pm, Km_func)
    realPart =  @. real(w_pm / sqrt_eval)        
    integrand = realPart .* evalKernel

    Zval = 1 .- trapz(w_prime, transpose(integrand)) ./ w_static

    # Evaluate Delta
    eval_real = @. real(Delta_func_eval / sqrt_eval)

    evalKernel = evaluate_Kernels(w_static, w_pm, Kp_func)
    int1 = @. eval_real * evalKernel
    int2 = @. mu *eval_real * tanh(beta * w_pm / 2)

    Delta = (trapz(w_prime, transpose(int1)) .- trapz(w_prime, transpose(int2))) ./ Zval

    return Zval, Delta, shift
    
end



"""
     realEliashbergEq(mu, beta, Delta_func_eval, sqrt_eval, Kp_func, Km_func, w_prime, w_max, W_pm, w_static, W_sm, W_prime) 

Real axis Eliashberg equations in cDOS+μ approximation
"""
function realEliashbergEq(mu::Float64, beta::Float64, Delta_func_eval::Vector{ComplexF64}, sqrt_eval::Vector{ComplexF64}, 
                          Kp_func::AbstractInterpolation, Km_func::AbstractInterpolation, w_prime::Vector{Float64}, 
                          w_static::StepRangeLen, w_pm::Vector{Float64}) 

    # Evaluate Z
    evalKernel = evaluate_Kernels(w_static, w_pm, Km_func)
    realPart =  @. real(w_pm / sqrt_eval)        
    integrand = @. realPart * evalKernel

    Zval = 1 .- trapz(w_prime, transpose(integrand)) ./ w_static

    # Evaluate Delta
    eval_real = @. real(Delta_func_eval / sqrt_eval)

    evalKernel = evaluate_Kernels(w_static, w_pm, Kp_func)
    int1 = @. eval_real * evalKernel
    int2 = @. mu *eval_real * tanh(beta * w_pm / 2)

    Delta = (trapz(w_prime, transpose(int1)) .- trapz(w_prime, transpose(int2))) ./ Zval

    return Zval, Delta

end


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
