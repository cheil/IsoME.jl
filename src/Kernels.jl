"""
    File containing the functions for the real axis solver
"""

function make_integration_axis(W_width, pts_cheb, pts_lin, epsilon)

    # Chebyshev grid inside [-epsilon, epsilon] to handle 1/Omega singularities.
    cheb_idx = div(pts_cheb, 2):(pts_cheb - 1)
    θ = @. (2 * cheb_idx + 1) / (2 * pts_cheb) * π

    W_cheb_right = @. epsilon * (1 + cos(θ))
    W_cheb_right = reverse(W_cheb_right)
    W_cheb_left = @. -epsilon * (1 + cos(θ))

    # Linear grids outside [-epsilon, epsilon].
    W_lin_left = range(-W_width, stop=-epsilon, length=div(pts_lin, 2))
    W_lin_right = range(epsilon, stop=W_width, length=div(pts_lin, 2))

    return collect(vcat(W_lin_left, W_cheb_left, W_cheb_right, W_lin_right))
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
    # w_axis must be same as w_static_chi
    if inp.cDOS_flag == 1
        w_axis = range(0, stop=inp.reOmega_c, length=5000)
    else
         # ω linear grid to store K(ω,ω'), has to be dimensions max(reOmega_c, reOmega_c_shift) + same for num
        w_axis = range(0, stop=max(inp.reOmega_c, inp.reOmega_c_shift), length=max(inp.numReal_c_shift,inp.numReal_c))

        # TEST
        #w_axis = range(0, stop=inp.reOmega_c, length=5000)
    end
    int_axis = make_integration_axis(W_right-W_left, 300, 300, 3)      # Ω-integration axis


    return W_cut, W_left, w_axis, int_axis
end

function kernel_integral_helper(
    s_rel::StepRangeLen{Float64},
    int_axis::Vector{Float64},
    W_cut::Float64,
    W_left::Float64,
    G::AbstractInterpolation,
    n,
    log_file=nothing;
    progress_counter=nothing,
    progress_interval=nothing,
    progress_state=nothing,
)

    N = length(s_rel)
    M = length(int_axis)
    integrand_1 = zeros(N, M)
    integrand_2 = zeros(N, M)
    @inbounds for i in 1:N
        for j in 1:M
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

        if !isnothing(progress_counter) && !isnothing(progress_interval)
            progress_counter[] += 1
            if mod(progress_counter[], progress_interval) == 0
                print_kernel_progress(log_file, progress_state)
            end
        end
    end

    integrand_1[isnan.(integrand_1)] .= 0.0
    integrand_2[isnan.(integrand_2)] .= 0.0

    integral_1 = trapz(int_axis, integrand_1)
    integral_2 = trapz(int_axis, integrand_2)

    return integral_1, integral_2
end


function precompute_integrals(num_w, w_max, W_left, W_cut, int_axis, G, n, log_file=nothing; progress_state=nothing)
    W_cut = W_cut[2]-W_cut[1]
    dw = w_max / (num_w - 1)
    progress_counter = isnothing(log_file) || isnothing(progress_state) ? nothing : Ref(0)
    progress_interval = isnothing(progress_counter) ? nothing : kernel_progress_interval(3 * (2 * num_w - 1), progress_state)

    s1_idx = 0:2*(num_w-1)
    s1_rel = (s1_idx .- 2*(num_w-1)) .* dw .- W_left
    I1_1, I1_2 = kernel_integral_helper(
        s1_rel,
        int_axis,
        W_cut,
        W_left,
        G,
        n,
        log_file;
        progress_counter=progress_counter,
        progress_interval=progress_interval,
        progress_state=progress_state,
    )

    s2_idx = 0:2*(num_w-1)
    s2_rel = (s2_idx .- (num_w-1)) .* dw .- W_left
    I2_1, I2_2 = kernel_integral_helper(
        s2_rel,
        int_axis,
        W_cut,
        W_left,
        G,
        n,
        log_file;
        progress_counter=progress_counter,
        progress_interval=progress_interval,
        progress_state=progress_state,
    )

    s4_idx = 0:2*(num_w-1)
    s4_rel = s4_idx .* dw .- W_left
    I4_1, I4_2 = kernel_integral_helper(
        s4_rel,
        int_axis,
        W_cut,
        W_left,
        G,
        n,
        log_file;
        progress_counter=progress_counter,
        progress_interval=progress_interval,
        progress_state=progress_state,
    )

    return [[I1_1, I1_2] [I2_1, I2_2] [I2_1, I2_2] [I4_1, I4_2]]
end

function compute_real_kernels(w_axis, num_w, w_max, W_left, W_cut, int_axis, f, n, G, log_file=nothing, line_width=80)

    println("num_w: ",num_w)
    println("w_max: ", w_max)

    progress_state = nothing
    if !isnothing(log_file)
        progress_state = start_kernel_progress("real", log_file, line_width)
    end

    w2_grid = repeat(w_axis, 1, length(w_axis))
    i1_grid = repeat(transpose((1:num_w)), num_w)
    i2_grid = repeat(1:num_w, 1, num_w)
    integrals = precompute_integrals(num_w, w_max, W_left, W_cut, int_axis, G, n, log_file; progress_state=progress_state) #slow

    if !isnothing(log_file)
        finish_kernel_progress(log_file, progress_state)
    end

    idx_a = @. -i1_grid - i2_grid + 2*num_w+1       #-w-w'
    idx_b = @. i1_grid - i2_grid + num_w            #w-w'
    idx_c = @. i2_grid - i1_grid + num_w            #w'-w
    idx_d = @. i1_grid + i2_grid - 1                #w+w'

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

function compute_imag_kernels(w_axis, f, n, a2F_itp, log_file=nothing, line_width=80)

    N = length(w_axis)

    progress_state = nothing
    if !isnothing(log_file)
        progress_state = start_kernel_progress("imag", log_file, line_width)
    end
    progress_interval = isnothing(log_file) ? nothing : kernel_progress_interval(N, progress_state)

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
            im_part2[i, j] = a2F_itp(dw) * (f(-w2) + n(dw))
        end

        if dw < 0  
            im_part3[i, j] = -a2F_itp(-dw) * (f_w2 + n(-dw))
        end
 
        im_part4[i, j] = a2F_itp(pw) * (f_w2 + n(pw))

        if !isnothing(log_file) && mod(i, progress_interval) == 0 && j == N
            print_kernel_progress(log_file, progress_state)
        end
    end

    Kp_imag = π .* (im_part2 .- im_part3 .- im_part4)
    Km_imag = π .* (-im_part2 .+ im_part3 .- im_part4)

    if !isnothing(log_file)
        finish_kernel_progress(log_file, progress_state)
    end

    return Kp_imag, Km_imag
end

# """
#     kernel_integral_helper(s_rel, int_axis, W_cut, W_left, a2F_itp)

# Compute the integrals I1(s_rel) and I2(s_rel) for the relative frequency coordinate s_rel = +- (w-w')
# int_axis is the integration axis with a denser grid around the singularity at x=0                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
# W_cut contains the upper and lower cutoffs of the Eliashberg spectral function a2F_itp
# """
# function kernel_integral_helper(s_rel::StepRangeLen{Float64}, int_axis::Vector{Float64}, W_cut::Vector{Float64}, a2F_itp::AbstractInterpolation, n)

#     W_width = W_cut[2] - W_cut[1]
#     N = length(s_rel)
#     M = length(int_axis)
#     integrand_1 = zeros(N, M)
#     integrand_2 = zeros(N, M)
#     @inbounds for i in 1:N, j in 1:M
#         s = s_rel[i]
#         wz = int_axis[j]                      
#         included = s > 0 && s < W_width

#         if included
#             x = wz + s + W_cut[1]
#             G = a2F_itp(x)
#             if G > 0
#                 integrand_2[i, j] = G / wz
#                 integrand_1[i, j] = G * n(x) / wz
#             end
#         else
#             x = wz + W_cut[1]
#             G = a2F_itp(x)
#             if G > 0
#                 integrand_2[i, j] = G / (wz - s)
#                 integrand_1[i, j] = G * n(x) / (wz - s)
#             end
#         end
#     end

#     integrand_1[isnan.(integrand_1)] .= 0.0
#     integrand_2[isnan.(integrand_2)] .= 0.0

#     integral_1 = trapz(int_axis, integrand_1)
#     integral_2 = trapz(int_axis, integrand_2)


#     return integral_1, integral_2
# end


# function precompute_integrals(num_w, w_axis, W_cut, int_axis, a2F_itp, n, W_left)
#     dw = w_axis[2] - w_axis[1]

#     println("cut", W_cut)
#     println("left", W_left)

#     s_idx = (-2*num_w-1):(2*num_w-1) 
#     s_rel = s_idx  .* dw .-W_left
#     I1, I2 = kernel_integral_helper(s_rel, int_axis, W_cut, a2F_itp, n)

#     return I1, I2
# end

# """
#     compute_real_kernels(w_axis, num_w, w_max, W_cut, int_axis, f, n, a2F_itp)

# Compute the real part of the kernel K(w,w') according to https://doi.org/10.1103/PhysRevB.54.6648
# Kp (Km) occurs in the integration for Δ (Z)
# """
# function compute_real_kernels(w_axis, num_w, w_max, W_cut, int_axis, f, n, a2F_itp, W_left)

#     num_w = length(w_axis)
#     println(num_w)
#     w_max = maximum(w_axis)

#     w2_grid = repeat(w_axis, 1, length(w_axis))
#     i1_grid = repeat((0:num_w-1)', num_w)
#     i2_grid = repeat(0:num_w-1, 1, num_w)
#     s_offset = 2*(num_w-1) +1
#     I1, I2 = precompute_integrals(num_w, w_axis, W_cut, int_axis, a2F_itp, n, W_left) 

#     Evalf_minus = f.(-w2_grid)
#     Evalf_plus = f.(w2_grid)

#     idx_a = @. -i1_grid - i2_grid + s_offset
#     idx_b = @. i1_grid - i2_grid + s_offset
#     idx_c = @. i2_grid - i1_grid + s_offset
#     idx_d = @. i1_grid + i2_grid + s_offset
#     re1 = transpose(Evalf_minus .* I2[idx_a] .+ I1[idx_a]) 
#     re2 = transpose(Evalf_minus .* I2[idx_b] .+ I1[idx_b])
#     re3 = transpose(Evalf_plus .* I2[idx_c] .+ I1[idx_c])
#     re4 = transpose(Evalf_plus .* I2[idx_d] .+ I1[idx_d])

#     Kp_real = re1 .+ re2 .- re3 .- re4
#     Km_real = re1 .- re2 .+ re3 .- re4

#     println("Km ", Km_real[1,1])

#     return Kp_real, Km_real
# end



"""
    kernels(β, inp, w_axis, W_cut, int_axis, a2F_itp)

Compute the kernels at a given temperature
Kp(ω,ω') = -K(ω,ω') + K(ω,-ω')
Km(ω,ω') = K(ω,ω') + K(ω,-ω')
"""
function kernels(β, inp, w_axis, W_cut, int_axis, a2F_itp, W_left, log_file=nothing, line_width=80)

    # Fermi-Dirac and Bose-Einstein distributions
    f = x -> 1 / (exp(β * x) + 1)
    n = x -> x == 0 ? 0.0 : abs(1 / (exp(β * x) - 1))

    # Compute real and imaginary parts of kernels
    Kp_imag, Km_imag = compute_imag_kernels(w_axis, f, n, a2F_itp, log_file, line_width)
    Kp_real, Km_real = compute_real_kernels(w_axis, length(w_axis), w_axis[end], W_left, W_cut, int_axis, f, n, a2F_itp, log_file, line_width)     # reOmega_c, numReal_c

    # Combine into complex kernels
    Kp = Kp_real .+ im .* (Kp_imag)
    Km = Km_real .+ im .* (Km_imag)

    # Kp_func = extrapolate(scale(interpolate(Kp, BSpline(Linear())), w_axis, w_axis), Flat())
    # Km_func = extrapolate(scale(interpolate(Km, BSpline(Linear())), w_axis, w_axis), Flat())

    Kp_func = interpolate((w_axis, w_axis), Kp, Gridded(Linear()))
    Km_func = interpolate((w_axis, w_axis), Km, Gridded(Linear()))

    return Kp_func, Km_func
end

function precompute(β, inp, matval, a2F_itp, console, log_file)

    printTextCentered("Precomputing Kernels", console["cDOS"]["partingLine"], file=log_file, bold=true)
    print_kernel_log(log_file, "\n")

    ### Set up axis and parameters ###
    W_cut, W_left, w_axis, int_axis = setUpAxis(inp, matval)

    Kp_func, Km_func = kernels(β, inp, w_axis, W_cut, int_axis, a2F_itp, W_left, log_file, length(console["cDOS"]["partingLine"]))

    w_static = range(1e-1, stop=inp.reOmega_c, length=inp.numReal_c)        # grid of Z(w), Delta(w), sensitive to start value, do not chose < 1e-1
    w_static_chi = range(1e-1, stop=inp.reOmega_c_shift, length=inp.numReal_c_shift)   # grid of χ(ω), 10*inp.reOmega_c
    w_dynam = reverse(inp.reOmega_c .+ inp.reOmega_c .* cos.((2 .* (inp.n_cheb/2:inp.n_cheb-1) .+ 1) .* π ./ (2 * inp.n_cheb)))    # only positvie chebyshev nodes

    # if inp.plot_flag
    #     plot_kernel_diagonal(β, inp, Kp_func, Km_func, w_static_chi, log_file)
    # end

    return (Kp_func, Km_func, w_static, w_dynam, w_static_chi)

end


###########################################
# --------------- Helpers --------------- #
###########################################
function print_kernel_log(log_file, text)
    if isnothing(log_file)
        print(text)
    else
        printTee(log_file, text)
    end
end

function start_kernel_progress(label, log_file, line_width)
    prefix = label * " Kernel: "
    progress_state = (
        printed=Ref(0),
        target=max(1, line_width - length(prefix) - 3),
        spinner=Ref(0),
        prefix=prefix,
        log_file=log_file,
    )
    redraw_kernel_progress(progress_state)
    return progress_state
end

function kernel_progress_interval(total_steps, progress_state)
    return max(1, cld(total_steps, progress_state.target))
end

function redraw_kernel_progress(progress_state)
    spinner = raw"-\|/"[mod(progress_state.spinner[], 4) + 1]
    progress_state.spinner[] += 1
    done = progress_state.printed[]
    remaining = progress_state.target - done
    print("\r" * progress_state.prefix * string(spinner) * "" *  "="^done * " "^remaining * "|")
end

function print_kernel_progress(log_file, progress_state)
    if progress_state.printed[] < progress_state.target
        progress_state.printed[] += 1
        redraw_kernel_progress(progress_state)
    end
end

function finish_kernel_progress(log_file, progress_state)
    progress_state.printed[] = progress_state.target
    progress_state.spinner[] = 2
    redraw_kernel_progress(progress_state)
    print("\n")
    if !isnothing(log_file)
        print(log_file, progress_state.prefix * "|" *  "="^progress_state.target * "|\n")
    end
end

function plot_kernel_diagonal(β, inp, Kp_func, Km_func, w_static_chi, log_file)
    try
        w_diag = collect(w_static_chi)
        Kp_diag = [Kp_func(w, w) for w in w_diag]
        Km_diag = [Km_func(w, w) for w in w_diag]

        itemp = round(1 / (kb * β), digits=4)
        temp_label = replace(string(itemp), "." => "p")

        plot(
            w_diag,
            real.(Kp_diag),
            label="Re K+",
            xlabel="ω / meV",
            ylabel="K(ω,ω)",
            linewidth=2,
        )
        plot!(w_diag, imag.(Kp_diag), label="Im K+", linewidth=2)
        plot!(w_diag, real.(Km_diag), label="Re K-", linewidth=2)
        plot!(w_diag, imag.(Km_diag), label="Im K-", linewidth=2)

        savefig("kernel_diagonal_T$(temp_label)K.png")
    catch ex
        printWarning("Error while plotting real-axis kernel diagonal.", log_file, ex=ex)
    end
end
