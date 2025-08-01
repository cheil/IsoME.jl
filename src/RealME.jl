"""
    File containing the functions for the real axis solver
"""

function make_kernel_grid(w_max, pts_lin)
    return range(0, stop=w_max, length=pts_lin)
end

function make_integration_axis(W_left, W_right, pts_cheb, pts_lin, epsilon)
    W_cut = (W_right - W_left)

    θ = (2*(div(pts_cheb,2):(pts_cheb-1)) .+ 1) .* (π / (2 * pts_cheb))
    W_cheb_right = epsilon .* (1 .+ cos.(θ)) |> reverse
    W_cheb_left = -epsilon .* (1 .+ cos.(θ))
    W_lin_left = range(-W_cut, stop=-epsilon, length=div(pts_lin,2))
    W_lin_right = range(epsilon, stop=W_cut, length=div(pts_lin,2))
    return vcat(W_lin_left, W_cheb_left, W_cheb_right, W_lin_right)
end

function setUpAxis(inp, matval)
    (a2f_omega, a2f) = matval

    # G(Ω): Interpolated a2f function
    G = linear_interpolation(a2f_omega, a2f, extrapolation_bc=Flat())

    # nonzero values of a2F, ensure endpoints are zero
    idx_left = findfirst(a2f .> 0) -1   
    idx_right = findlast(a2f .> 0) +1

    W_left = a2f_omega[idx_left]
    W_right = a2f_omega[idx_right]

    # Parameters 
    W_cut = W_right - W_left
    # Frequency and integration grids
    w_axis = make_kernel_grid(inp.real_c, inp.numReal_c)
    int_axis = make_integration_axis(W_left, W_right, 300, 300, 3)

    return W_cut, w_axis, int_axis, W_left, G
end

function kernel_integral_helper(s_rel::StepRangeLen{Float64}, int_axis::Vector{Float64}, W_cut::Float64, W_left::Float64, G::AbstractInterpolation, n)

    N = length(s_rel)
    M = length(int_axis)
    integrand_1 = zeros(N, M)
    integrand_2 = zeros(N, M)
    @inbounds for i in 1:N, j in 1:M
        s = s_rel[i]
        wz = int_axis[j]

        included = s > 0 && s < W_cut
        if included
            G1 = G(wz + s + W_left)
            if G1 > 0
                integrand_1[i, j] = G1 / wz

                integrand_2[i, j] = integrand_1[i, j] * n(wz + s + W_left)
            end
        else
            G2 = G(wz + W_left)
            if G2 > 0 
                integrand_1[i, j] = G2 / (wz - s)
 
                integrand_2[i, j] = integrand_1[i, j] * n(wz + W_left)
            end
        end
    end

    integrand_1[isnan.(integrand_1)] .= 0.0
    integrand_2[isnan.(integrand_2)] .= 0.0

    integral_1 = trapz(int_axis, integrand_1)
    integral_2 = trapz(int_axis, integrand_2)


    return integral_1, integral_2
end


function precompute_integrals(num_w, w_max, W_left, W_cut, int_axis, G, n)
    dw = w_max / (num_w - 1)

    s1_idx = 0:2*(num_w-1)
    s1_rel = (s1_idx .- 2*(num_w-1)) .* dw .- W_left
    I1_1, I1_2 = kernel_integral_helper(s1_rel, int_axis, W_cut, W_left, G, n)

    s2_idx = 0:2*(num_w-1)
    s2_rel = (s2_idx .- (num_w-1)) .* dw .- W_left
    I2_1, I2_2 = kernel_integral_helper(s2_rel, int_axis, W_cut, W_left, G, n)

    s4_idx = 0:2*(num_w-1)
    s4_rel = s4_idx .* dw .- W_left
    I4_1, I4_2 = kernel_integral_helper(s4_rel, int_axis, W_cut, W_left, G, n)

    return [[I1_1, I1_2] [I2_1, I2_2] [I2_1, I2_2] [I4_1, I4_2]]
end

function compute_real_kernels(w_axis, num_w, w_max, W_left, W_cut, int_axis, f, n, G)

    w2_grid = repeat(w_axis, 1, length(w_axis))
    i1_grid = repeat((1:num_w)', num_w)
    i2_grid = repeat(1:num_w, 1, num_w)
    integrals = precompute_integrals(num_w, w_max, W_left, W_cut, int_axis, G, n) #slow

    idx_a = @. -i1_grid - i2_grid + 2*num_w+1
    idx_b = @. i1_grid - i2_grid + num_w
    idx_c = @. i2_grid - i1_grid + num_w
    idx_d = @. i1_grid + i2_grid - 1

    Evalf_minus = f.(-w2_grid)
    Evalf_plus = f.(w2_grid)
    re1 = transpose(Evalf_minus .* integrals[1, 1][idx_a] .+ integrals[2, 1][idx_a]) 
    re2 = transpose(Evalf_minus .* integrals[1, 2][idx_b] .+ integrals[2, 2][idx_b])
    re3 = transpose(Evalf_plus .* integrals[1, 3][idx_c] .+ integrals[2, 3][idx_c])
    re4 = transpose(Evalf_plus .* integrals[1, 4][idx_d] .+ integrals[2, 4][idx_d])

    Kp_real = re1 .+ re2 .- re3 .- re4
    Km_real = re1 .- re2 .+ re3 .- re4

    return Kp_real, Km_real
end

function compute_imag_kernels(w_axis, f, n, G)

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
            im_part2[i, j] = G(dw) * (1-f_w2 + n(dw))
        end

        if -dw > 0  
            im_part3[i, j] = G(-dw) * (f_w2 + n(-dw))
        end
 
        im_part4[i, j] = G(pw) * (f_w2 + n(pw))
    end

    Kp_imag = π .* (im_part2 .- im_part3 .- im_part4)
    Km_imag = π .* (-im_part2 .+ im_part3 .- im_part4)

    return Kp_imag, Km_imag
end

"""
    kernels()

Compute kernels at given temperature
"""
function kernels(β, inp, w_axis, W_left, W_cut, int_axis, G)

    # Fermi-Dirac and Bose-Einstein distributions
    f = x -> 1 / (exp(β * x) + 1)
    n = x -> x == 0 ? 0.0 : abs(1 / (exp(β * x) - 1))

    # Compute real and imaginary parts of kernels
    Kp_imag, Km_imag = compute_imag_kernels(w_axis, f, n, G)
    Kp_real, Km_real = compute_real_kernels(w_axis, inp.numReal_c, inp.real_c, W_left, W_cut, int_axis, f, n, G)

    # Combine into complex kernels
    Kp = Kp_real .+ im .* Kp_imag
    Km = Km_real .+ im .* Km_imag

    Kp_func = interpolate((w_axis, w_axis), Kp, Gridded(Linear()))
    Km_func = interpolate((w_axis, w_axis), Km, Gridded(Linear()))


    return Kp_func, Km_func
end

function precompute(β, inp, w_axis, W_left, W_cut, int_axis, G)

    Kp_func, Km_func = kernels(β, inp, w_axis, W_left, W_cut, int_axis, G)

     # new cheb grid
    w_static = range(1e-4, stop=inp.real_c, length=inp.numReal_c)
    w_dynam = reverse(inp.real_c .+ inp.real_c .* cos.((2 .* ((inp.n_cheb/2):(inp.n_cheb-1)) .+ 1) .* π ./ (2 * inp.n_cheb)))
    W_static = repeat(transpose(w_static), length(w_dynam))
    W_dynam = repeat(w_dynam, 1, length(w_static))

    return (Kp_func, Km_func, w_static, w_dynam, W_static, W_dynam)

end



