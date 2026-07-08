"""
    File containing the functions for the real axis solver
"""

"""
    make_integration_axis(W_width, pts_cheb, pts_lin, epsilon)

Set up the Ω-integration axis
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
    setUpOmegaAxis(inp, mataval)

Set up the integration axis for the Ω integration
"""
function setUpOmegaAxis(inp, matval)
    (a2f_omega, a2f) = matval


    # nonzero values of a2F, ensure endpoints are zero
    idx_left = findfirst(a2f .> 1e-6)   
    idx_right = findlast(a2f .> 1e-6)

    W_left = a2f_omega[idx_left]
    W_right = a2f_omega[idx_right]

    # interpolate
    G = scale(interpolate(a2f[idx_left:idx_right], BSpline(Linear())), a2f_omega[idx_left:idx_right])

    # Frequency and integration grids
    # w_axis must be same as w_static_chi
    if inp.cDOS_flag == 1
        w_axis = 0:inp.dOmega:inp.reOmega_c   # w-axis at which K(ω,ω') is calculated
        #w_axis = range(0, stop=inp.reOmega_c, length=5000) 
    else
        # ω linear grid to store K(ω,ω'), has to span max(reOmega_c, reOmega_c_shift).
        # w_axis must match w_static_chi, so it uses the χ step size domega_shift.
        w_axis = 0:inp.dOmega:max(inp.reOmega_c, inp.reOmega_c_shift)
    end
    int_axis = make_integration_axis(W_right-W_left, 300, 300, 3)      # Ω-integration axis #


    return G, W_left, W_right, w_axis, int_axis
end

"""
    pv_linear(g, x, E)

Cauchy principal-value integral  P ∫ g(Ω)/(Ω − E) dΩ  of a piecewise-linear g
(knots `x`, values `g`), evaluated in closed form segment by segment (Method B):

    ∫_{Ω_k}^{Ω_{k+1}} (aΩ+b)/(Ω-E) dΩ = a·(Ω_{k+1}-Ω_k) + (aE+b)·ln|(Ω_{k+1}-E)/(Ω_k-E)|.

Exact for piecewise-linear g, valid for E inside or outside [x[1], x[end]], and needs no
integration grid or symmetric-cancellation trick. A ln-argument that hits a knot exactly is
dropped: the coefficient (aE+b) is continuous across the shared knot, so that singular
contribution → 0 in the limit.
"""
function pv_linear(g::Vector{Float64}, x::Vector{Float64}, E::Float64)
    acc = 0.0
    @inbounds for k in 1:length(x)-1
        dx = x[k+1] - x[k]
        a  = (g[k+1] - g[k]) / dx
        b  = g[k] - a * x[k]
        d1 = x[k+1] - E
        d2 = x[k]   - E
        l1 = abs(d1) < 1e-12 ? 0.0 : log(abs(d1))
        l2 = abs(d2) < 1e-12 ? 0.0 : log(abs(d2))
        acc += a * dx + (a * E + b) * (l1 - l2)
    end
    return acc
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
                if wz+s > 0 && wz+s < W_cut # G = 0 outside
                    G1 = G(wz + s + W_left)
                    integrand_1[i, j] = G1 / wz
                    integrand_2[i, j] = integrand_1[i, j] * n(wz + s + W_left)
                end
            else

                if wz > 0 && wz < W_cut   # G = 0 outside
                    G2 = G(wz + W_left)
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

    println("nan: ", sum(isnan.(integrand_1)))
    integrand_1[isnan.(integrand_1)] .= 0.0
    integrand_2[isnan.(integrand_2)] .= 0.0


    integral_1 = trapz(int_axis, integrand_1)
    integral_2 = trapz(int_axis, integrand_2)

   
    # ---- Method B: exact analytic PV of the piecewise-linear α²F, for comparison ----
    # Same integrals  I₁(s)=P∫α²F(Ω)/(Ω-E)dΩ  and  I₂(s)=P∫α²F(Ω)n(Ω)/(Ω-E)dΩ , E=W_left+s,
    # but done in closed form over the band [W_left, W_right] instead of the trapz-over-int_axis
    # scheme above. Reports the largest A(trapz) − B(analytic) discrepancy over the s-vector.
    
    # print("B: ")
    # @time begin
    # Ω_grid = collect(range(W_left, W_right, length=2000))
    # g1 = G.(Ω_grid); g1[g1 .< 0.0] .= 0.0
    # g2 = g1 .* n.(Ω_grid)
    # integral_1_B = zeros(N)
    # integral_2_B = zeros(N)
    # @inbounds for i in 1:N
    #     E = W_left + s_rel[i]
    #     integral_1_B[i] = pv_linear(g1, Ω_grid, E)
    #     integral_2_B[i] = pv_linear(g2, Ω_grid, E)
    # end
    # d1v = abs.(integral_1 .- integral_1_B)
    # d2v = abs.(integral_2 .- integral_2_B)
    # i1m = argmax(d1v); i2m = argmax(d2v)
    # @info "kernel PV  A(trapz) − B(analytic)" maxΔI1=d1v[i1m] at_s1=s_rel[i1m] maxΔI2=d2v[i2m] at_s2=s_rel[i2m]
    # end
    # ------------------------------------------------------------------------------------------

    return integral_1, integral_2
end


function precompute_integrals(num_w, w_max, W_left, W_right,int_axis, G, n, log_file=nothing; progress_state=nothing)
    W_cut = W_right - W_left
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

function compute_real_kernels(w_axis, num_w, w_max, W_left, W_right, int_axis, f, n, G, log_file=nothing, line_width=80)

    progress_state = nothing
    if !isnothing(log_file)
        progress_state = start_kernel_progress(" Re", log_file, line_width)
    end

    w2_grid = repeat(w_axis, 1, length(w_axis))
    i1_grid = repeat(transpose((1:num_w)), num_w)
    i2_grid = repeat(1:num_w, 1, num_w)
    integrals = precompute_integrals(num_w, w_max, W_left, W_right, int_axis, G, n, log_file; progress_state=progress_state) #slow

    println("Hallo: ", maximum(maximum(integrals)))

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

    println("Km real: ", maximum(abs.(re1)))
    println("Km real: ", maximum(abs.(re2)))
    println("Km real: ", maximum(abs.(re3)))
    println("Km real: ", maximum(abs.(re4)))
    
    println("Km real: ", maximum(Km_real))

    return Kp_real, Km_real
end


function compute_imag_kernels(w_axis, f, n, G, W_left, W_right, log_file=nothing, line_width=80)


    N = length(w_axis)

    progress_state = nothing
    if !isnothing(log_file)
        progress_state = start_kernel_progress(" Im", log_file, line_width)
    end
    progress_interval = isnothing(log_file) ? nothing : kernel_progress_interval(N, progress_state)

    Kp_imag = zeros(N,N)
    Km_imag = zeros(N,N)
    @inbounds for j in 1:N
        w2 = w_axis[j]     
        for i in 1:N
            w1 = w_axis[i]  
  
            dw = w1 - w2
            pw = w1 + w2
            # G = 0 outside [W_left, W_right]
            if dw > W_left && dw < W_right
                temp = G(dw) * (f(-w2)+ n(dw))
                Kp_imag[i, j] += temp
                Km_imag[i, j] -= temp
            elseif dw > -W_right && dw < -W_left
                temp = -G(-dw) * (f(w2) + n(-dw))
                Kp_imag[i, j] -= temp
                Km_imag[i, j] += temp
            end
    
            if pw > W_left && pw < W_right
                temp = G(pw) * (f(w2) + n(pw))
                Kp_imag[i, j] -= temp
                Km_imag[i, j] -= temp
            end

            if !isnothing(log_file) && mod(i, progress_interval) == 0 && j == N
                print_kernel_progress(log_file, progress_state)
            end
        end
    end


    @. Kp_imag *= π
    @. Km_imag *= π


    return Kp_imag, Km_imag
end


"""
    kernels(β, inp, w_axis, W_cut, int_axis, G)

Compute the kernels at a given temperature
Kp(ω,ω') = -K(ω,ω') + K(ω,-ω')
Km(ω,ω') = K(ω,ω') + K(ω,-ω')
"""
function kernels(β, inp, w_axis, W_left, W_right, int_axis, G, log_file=nothing, line_width=80)

    # Fermi-Dirac and Bose-Einstein distributions
    f = x -> 1 / (exp(β * x) + 1)
    n = x -> x == 0 ? 0.0 : abs(1 / (exp(β * x) - 1))

    # Compute real and imaginary parts of kernels
    @time Kp_imag, Km_imag = compute_imag_kernels(w_axis, f, n, G, W_left, W_right, log_file, line_width)
    @time Kp_real, Km_real = compute_real_kernels(w_axis, length(w_axis), w_axis[end], W_left, W_right, int_axis, f, n, G, log_file, line_width)     # reOmega_c, domega

    # Combine into complex kernels
    Kp = Kp_real .+ im .* (Kp_imag)
    Km = Km_real .+ im .* (Km_imag)

    # more accurate but slower
    # Kp_func = interpolate((w_axis, w_axis), Kp, Gridded(Linear()))
    # Km_func = interpolate((w_axis, w_axis), Km, Gridded(Linear()))

    # faster
    Kp_func = scale(interpolate(Kp, BSpline(Linear())), w_axis, w_axis)
    Km_func = scale(interpolate(Km, BSpline(Linear())), w_axis, w_axis)

    return Kp_func, Km_func
end

function precompute(β, inp, matval, console, log_file)

    printTextCentered("Precomputing Kernels", console["cDOS"]["partingLine"], file=log_file, bold=true)
    print_kernel_log(log_file, "\n")

    ### Set up axis and parameters ###
    G, W_left, W_right, w_axis, int_axis = setUpOmegaAxis(inp, matval)

    Kp_func, Km_func = kernels(β, inp, w_axis, W_left, W_right, int_axis, G, log_file, length(console["cDOS"]["partingLine"]))

    w_static = 1e-1:inp.domega:inp.reOmega_c        # grid of Z(w), Delta(w), sensitive to start value, do not chose < 1e-1
    w_static_chi = 1e-1:inp.domega_shift:inp.reOmega_c_shift   # grid of χ(ω), 10*inp.reOmega_c
    w_dynam = reverse(inp.reOmega_c .+ inp.reOmega_c .* cos.((2 .* (inp.n_cheb/2:inp.n_cheb-1) .+ 1) .* π ./ (2 * inp.n_cheb)))    # only positvie chebyshev nodes
    #append!(w_dynam, reverse(-w_dynam))   # symmetric chebyshev nodes
    # assert grid sizes
    @assert first(w_axis) <= first(w_static) && last(w_static) <= last(w_axis)
    if inp.include_Weep == 1
        @assert first(w_axis) <= first(w_static_chi) && last(w_static_chi) <= last(w_axis)
    end

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
    prefix = label * " K(w,w'): "
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

