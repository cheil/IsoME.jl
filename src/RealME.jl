"""
    File containing the functions for the real axis solver
"""

function make_integration_axis(W_right, pts_cheb, pts_lin, epsilon)

    θ = (2*(div(pts_cheb,2):(pts_cheb-1)) .+ 1) .* (π / (2 * pts_cheb))
    W_cheb_right = epsilon .* (1 .+ cos.(θ)) |> reverse
    W_cheb_left = -epsilon .* (1 .+ cos.(θ))
    W_lin_left = range(-W_right, stop=-epsilon, length=div(pts_lin,2))
    W_lin_right = range(epsilon, stop=W_right, length=div(pts_lin,2))
    return vcat(W_lin_left, W_cheb_left, W_cheb_right, W_lin_right)
end



"""
    setUpAxis(inp, mataval)

Set up the integration axis for the Ω integration
"""
function setUpAxis(inp, matval)
    (a2f_omega, a2f) = matval

    # nonzero values of a2F, ensure endpoints are zero
    idx_left = findfirst(a2f .> 1e-6) -1   
    idx_right = findlast(a2f .> 1e-6) +1

    W_left = a2f_omega[idx_left]
    W_right = a2f_omega[idx_right]

    # Parameters 
    W_cut = [W_left, W_right]
    
    # Frequency and integration grids
    w_axis = range(0, stop=inp.reOmega_c, length=inp.numReal_c) # omega axis
    int_axis = make_integration_axis(W_right, 300, 300, 3) # Omega axis

    return W_cut, w_axis, int_axis
end

"""
    kernel_integral_helper(s_rel, int_axis, W_cut, W_left, a2F_itp)

Compute the integrals I1(s_rel) and I2(s_rel) for the relative frequency coordinate s_rel = +- (w-w')
int_axis is the integration axis with a denser grid around the singularity at x=0                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
W_cut contains the upper and lower cutoffs of the Eliashberg spectral function a2F_itp
"""
function kernel_integral_helper(s_rel::StepRangeLen{Float64}, int_axis::Vector{Float64}, W_cut::Vector{Float64}, a2F_itp::AbstractInterpolation, n)

    N = length(s_rel)
    M = length(int_axis)
    integrand_1 = zeros(N, M)
    integrand_2 = zeros(N, M)
    @inbounds for i in 1:N, j in 1:M
        s = s_rel[i]
        wz = int_axis[j]                      
        included = s > 0 && s < W_cut[2]

        if included
            if wz + s > W_cut[1] && wz+s < W_cut[2]
                G1 = a2F_itp(wz + s)
                integrand_2[i, j] = G1 / wz

                integrand_1[i, j] = integrand_2[i, j] * n(wz + s)
            end
        else
            if wz > W_cut[1] && wz < W_cut[2]
                G2 = a2F_itp(wz)
                integrand_2[i, j] = G2 / (wz - s)
 
                integrand_1[i, j] = integrand_2[i, j] * n(wz)
            end
        end
    end

    integrand_1[isnan.(integrand_1)] .= 0.0
    integrand_2[isnan.(integrand_2)] .= 0.0

    integral_1 = trapz(int_axis, integrand_1)
    integral_2 = trapz(int_axis, integrand_2)


    return integral_1, integral_2
end


function precompute_integrals(num_w, w_max, W_cut, int_axis, a2F_itp, n)
    dw = w_max / (num_w - 1)

    s_idx = -2*num_w:2*num_w
    s_rel = s_idx  .* dw
    I1, I2 = kernel_integral_helper(s_rel, int_axis, W_cut, a2F_itp, n)

    return I1, I2
end

"""
    compute_real_kernels(w_axis, num_w, w_max, W_cut, int_axis, f, n, a2F_itp)

Compute the real part of the kernel K(w,w') according to https://doi.org/10.1103/PhysRevB.54.6648
Kp (Km) occurs in the integration for Δ (Z)
"""
function compute_real_kernels(w_axis, num_w, w_max, W_cut, int_axis, f, n, a2F_itp)

    w2_grid = repeat(w_axis, 1, length(w_axis))
    i1_grid = repeat((1:num_w)', num_w)
    i2_grid = repeat(1:num_w, 1, num_w)
    I1, I2 = precompute_integrals(num_w, w_max, W_cut, int_axis, a2F_itp, n) 

    Evalf_minus = f.(-w2_grid)
    Evalf_plus = f.(w2_grid)

    idx_a = @. -i1_grid - i2_grid + 2*num_w+3
    idx_b = @. i1_grid - i2_grid + 2*num_w+1
    idx_c = @. i2_grid - i1_grid + 2*num_w+1
    idx_d = @. i1_grid + i2_grid + 2*num_w-1
    re1 = transpose(Evalf_minus .* I2[idx_a] .+ I1[idx_a]) 
    re2 = transpose(Evalf_minus .* I2[idx_b] .+ I1[idx_b])
    re3 = transpose(Evalf_plus .* I2[idx_c] .+ I1[idx_c])
    re4 = transpose(Evalf_plus .* I2[idx_d] .+ I1[idx_d])

    Kp_real = re1 .+ re2 .- re3 .- re4
    Km_real = re1 .- re2 .+ re3 .- re4

    return Kp_real, Km_real
end

function compute_imag_kernels(w_axis, f, n, a2F_itp)

    N = length(w_axis)

    # im_part1 = always 0
    im_part2 = zeros(N, N)
    im_part3 = zeros(N, N)
    im_part4 = zeros(N, N)

    @inbounds for i in 1:N, j in 1:N
        w1 = w_axis[i]  
        w2 = w_axis[j]          

        dw = w1 - w2
        pw = w1 + w2
        f_w2 = f(w2)

        if dw > 0
            im_part2[i, j] = a2F_itp(dw) * (1-f_w2 + n(dw))
        end

        if dw < 0  
            im_part3[i, j] = a2F_itp(-dw) * (f_w2 + n(-dw))
        end
 
        im_part4[i, j] = a2F_itp(pw) * (f_w2 + n(pw))
    end

    Kp_imag = π .* (im_part2 .- im_part3 .- im_part4)
    Km_imag = π .* (-im_part2 .+ im_part3 .- im_part4)


    return Kp_imag, Km_imag
end

"""
    kernels(β, inp, w_axis, W_cut, int_axis, a2F_itp)

Compute the Ω kernels at a given temperature
"""
function kernels(β, inp, w_axis, W_cut, int_axis, a2F_itp)

    # Fermi-Dirac and Bose-Einstein distributions
    f = x -> 1 / (exp(β * x) + 1)
    n = x -> x == 0 ? 0.0 : abs(1 / (exp(β * x) - 1))

    # Compute real and imaginary parts of kernels
    Kp_imag, Km_imag = compute_imag_kernels(w_axis, f, n, a2F_itp)
    Kp_real, Km_real = compute_real_kernels(w_axis, inp.numReal_c, inp.reOmega_c, W_cut, int_axis, f, n, a2F_itp)

    # Combine into complex kernels
    Kp = Kp_real .+ im .* Kp_imag
    Km = Km_real .+ im .* Km_imag

    #Kp_func = interpolate((w_axis, w_axis), Kp, Gridded(Linear()))
    #Km_func = interpolate((w_axis, w_axis), Km, Gridded(Linear()))

    Kp_func = scale(interpolate(Kp, BSpline(Linear())), w_axis, w_axis)
    Km_func = scale(interpolate(Km, BSpline(Linear())), w_axis, w_axis)


    return Kp_func, Km_func
end

function precompute(β, inp, w_axis, W_cut, int_axis, a2F_itp)

    Kp_func, Km_func = kernels(β, inp, w_axis, W_cut, int_axis, a2F_itp)

    # new cheb grid
    w_static = range(1e-4, stop=inp.reOmega_c, length=inp.numReal_c)        # grid of Z(w), Delta(w),...
    w_dynam = reverse(inp.reOmega_c .+ inp.reOmega_c .* cos.((2 .* (1:inp.n_cheb) .+ 1) .* π ./ (2 * inp.n_cheb)))

    return (Kp_func, Km_func, w_static, w_dynam)

end



