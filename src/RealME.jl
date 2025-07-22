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

    # Plot a2f
    idx_left = findfirst(a2f .> 1e-12)
    idx_right = findlast(a2f .> 1e-12)
    plot(a2f_omega, G(a2f_omega), xlabel="Ω [meV]", ylabel="α2F(Ω)", label="", lw=2)
    vline!([a2f_omega[idx_left], a2f_omega[idx_right]], label="", color=:red, show=true)
    savefig("a2f.png")

    W_left = a2f_omega[idx_left]
    W_right = a2f_omega[idx_right]

    # Parameters 
    W_cut = W_right - W_left
    # Frequency and integration grids
    w_axis = make_kernel_grid(inp.real_c, inp.numReal_c)
    int_axis = make_integration_axis(W_left, W_right, 300, 300, 3)

    return W_cut, w_axis, int_axis, W_left
end

function kernel_integral_helper(s_rel, int_axis, W_cut, W_left, G, n, flag = false)
    S = repeat(transpose(s_rel), length(int_axis))
    Wz = repeat(int_axis, 1, length(s_rel))
    included = (0 .< S) .& (S .< W_cut)

    integrand_1 = zeros(size(Wz))
    # First term: where singularity lies in the peak width
    EvalG1 = G.(Wz .+ S .+ W_left)
    mask1 = included .& (EvalG1 .> 0)
    EvalG1[mask1]
    integrand_1[mask1] .= EvalG1[mask1] ./ Wz[mask1]
    # Second term: where singularity lies outside peak width
    EvalG2 = G.(Wz .+ W_left)
    mask2 = .!included .& (EvalG2 .> 0)
    integrand_1[mask2] .= EvalG2[mask2] ./ (Wz[mask2] .- S[mask2])
    # Replace NaNs with 0
    integrand_1[isnan.(integrand_1)] .= 0.0

    # Integrate along Wz direction (i.e. axis=1 in NumPy → dimension 2 in Julia)
    integral_1 = trapz(int_axis, transpose(integrand_1))
#    integral_1 = map(i -> trapz(int_axis, integrand_1[i, :]), 1:size(integrand_1, 1))

    # Allocate result array
    integrand_2 = zeros(size(Wz))
    # Compute conditions
    mask1 = included .& (EvalG1 .> 0)
    mask2 = .!included .& (EvalG2 .> 0)
    Evaln1 = n.(Wz[mask1] .+ S[mask1] .+ W_left) 
    Evaln2 =  n.(Wz[mask2] .+ W_left)
    # Apply those conditions
    integrand_2[mask1] .= EvalG1[mask1] .* Evaln1 ./ Wz[mask1]
    integrand_2[mask2] .= EvalG2[mask2] .* Evaln2 ./ (Wz[mask2] .- S[mask2])
    
    # Replace any NaNs (e.g. from 0/0 or Inf) with 0.0
    integrand_2[isnan.(integrand_2)] .= 0.0
    # Integrate along the Wz axis (i.e., axis=1 in NumPy → dims=2 in Julia)
    integral_2 = trapz(int_axis, transpose(integrand_2))

    return integral_1, integral_2
end

function precompute_integrals(num_w, w_max, W_left, W_cut, int_axis, G, n)
    dw = w_max / (num_w - 1)

    s1_idx = 0:2*(num_w-1)
    s1_rel = (s1_idx .- 2*(num_w-1)) .* dw .- W_left
    I1_1, I1_2 = kernel_integral_helper(s1_rel, int_axis, W_cut, W_left, G, n, true)

    s2_idx = 0:2*(num_w-1)
    s2_rel = (s2_idx .- (num_w-1)) .* dw .- W_left
    I2_1, I2_2 = kernel_integral_helper(s2_rel, int_axis, W_cut, W_left, G, n)

    s4_idx = 0:2*(num_w-1)
    s4_rel = s4_idx .* dw .- W_left
    I4_1, I4_2 = kernel_integral_helper(s4_rel, int_axis, W_cut, W_left, G, n)

    return [[I1_1, I1_2] [I2_1, I2_2] [I2_1, I2_2] [I4_1, I4_2]]
end

function compute_real_kernels(w_axis, num_w, w_max, W_left, W_cut, int_axis, f, n, freq, a2F)
    # G(Ω): Interpolated a2F function, clamped at max(freq)
    G = linear_interpolation(freq, a2F, extrapolation_bc=Flat())

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

function compute_imag_kernels(w_axis, f, n, freq, a2F)
    # G(Ω): Interpolated a2F function, clamped at max(freq)
    G = linear_interpolation(freq, a2F, extrapolation_bc=Flat())

    N = length(w_axis)
    w1_grid = repeat(transpose(w_axis), N)
    w2_grid = repeat(w_axis, 1, N)

    im_part1 = zeros(N, N)
    im_part2 = zeros(N, N)
    im_part3 = zeros(N, N)
    im_part4 = zeros(N, N)

    for i in 1:N, j in 1:N
        w1 = w1_grid[i, j]
        w2 = w2_grid[i, j]
        if -w1 - w2 > 0
            im_part1[i, j] = G(-w1 - w2)*(f(-w2) + n(-w1 - w2))
        end
        if w1 - w2 > 0
            im_part2[i, j] = G(w1 - w2)*(f(-w2) + n(w1 - w2))
        end
        if w2 - w1 > 0
            im_part3[i, j] = G(w2 - w1)*(f(w2) + n(w2 - w1))
        end
        if w1 + w2 > 0
            im_part4[i, j] = G(w1 + w2)*(f(w2) + n(w1 + w2))
        end
    end
    im_part1 = transpose(im_part1)
    im_part2 = transpose(im_part2)
    im_part3 = transpose(im_part3)
    im_part4 = transpose(im_part4)

    Kp_imag = π .* (im_part1 .+ im_part2 .- im_part3 .- im_part4)
    Km_imag = π .* (im_part1 .- im_part2 .+ im_part3 .- im_part4)
    return Kp_imag, Km_imag
end

"""
    kernels()

Compute kernels at given temperature
"""
function kernels(β, inp, matval, w_axis, W_left, W_cut, int_axis)
    (a2f_omega, a2f) = matval

    # Fermi-Dirac and Bose-Einstein distributions
    f = x -> 1 / (exp(β * x) + 1)
    n = x -> x == 0 ? 0.0 : abs(1 / (exp(β * x) - 1))

    # Compute real and imaginary parts of kernels
    Kp_imag, Km_imag = compute_imag_kernels(w_axis, f, n, a2f_omega, a2f)
    Kp_real, Km_real = compute_real_kernels(w_axis, inp.numReal_c, inp.real_c, W_left, W_cut, int_axis, f, n, a2f_omega, a2f)

    # Combine into complex kernels
    Kp = Kp_real .+ im .* Kp_imag
    Km = Km_real .+ im .* Km_imag

    Kp_func = interpolate((w_axis, w_axis), Kp, Gridded(Linear()))
    Km_func = interpolate((w_axis, w_axis), Km, Gridded(Linear()))


    return Kp_func, Km_func
end

function Z(sqrt_eval, Km_func, W_prime, w_prime, w_static, W_sm, W_pm, w_max)

    integrand =  @. real(W_pm / sqrt_eval) * Km_func.(W_sm, W_pm)
    integrand[W_prime .>= w_max] .= 0

    temp = 1 .- trapz(w_prime, transpose(integrand)) ./ w_static

    return temp
end

function newDelta(mu, beta, Delta_func_eval, sqrt_eval, Z_val, Kp_func, w_prime, w_max, W_pm, W_sm, W_prime)

    eval_real = @. real(Delta_func_eval / sqrt_eval)

    int1 = @. eval_real * Kp_func(W_sm, W_pm)
    int2 = @. mu *eval_real * tanh(beta * W_pm / 2)

    int1[W_prime .>= w_max] .= 0
    int2[W_prime .>= w_max] .= 0

    temp = trapz(w_prime, transpose(int1)) .- trapz(w_prime, transpose(int2))
    
    return temp ./ Z_val
end

