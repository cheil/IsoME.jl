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

    #println("nan: ", sum(isnan.(integrand_1)))
    # I think can be removed as no nan's should happen anymore
    integrand_1[isnan.(integrand_1)] .= 0.0
    integrand_2[isnan.(integrand_2)] .= 0.0


    integral_1 = trapz(int_axis, integrand_1)
    integral_2 = trapz(int_axis, integrand_2)

    return integral_1, integral_2
end


function precompute_integrals(num_w, w_max, W_left, W_right,int_axis, G, n, log_file=nothing; progress_state=nothing)
    W_cut = W_right - W_left
    dw = w_max / (num_w - 1)
    progress_counter = isnothing(log_file) || isnothing(progress_state) ? nothing : Ref(0)
    progress_interval = isnothing(progress_counter) ? nothing : kernel_progress_interval(3 * (2 * num_w - 1), progress_state)

    # single unified s-grid: position p ↔ frequency index m = p - (2*num_w - 1),
    # covering every m ∈ -2(num_w-1) … 2(num_w-1) that the ±w±w' combinations reach.
    s_idx = -2*(num_w-1):2*(num_w-1)
    s_rel = s_idx .* dw .- W_left
    I1, I2 = kernel_integral_helper(
        s_rel,
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


    return I1, I2
end

function compute_real_kernels(w_axis, num_w, w_max, W_left, W_right, int_axis, f, n, G, log_file=nothing, line_width=80)

    progress_state = nothing
    if !isnothing(log_file)
        progress_state = start_kernel_progress(" Re", log_file, line_width)
    end

    I1, I2 = precompute_integrals(num_w, w_max, W_left, W_right, int_axis, G, n, log_file; progress_state=progress_state) #slow

    if !isnothing(log_file)
        finish_kernel_progress(log_file, progress_state)
    end

    N = length(w_axis)
    Kp_real = Matrix{Float64}(undef, N, N)
    Km_real = Matrix{Float64}(undef, N, N)
    @inbounds for j in 1:N, i in 1:N              # i innermost (column-major)
        w  = w_axis[j]                            
        fm = f(-w); fp = f(w)
        # indices into the single unified grid: p(m) = m + 2*num_w - 1, with
        # m = -(i+j-2) [-w-w'],  i-j [w-w'],  j-i [w'-w],  i+j-2 [w+w']
        ia = -i-j+2num_w+1; ib = i-j+2num_w-1; ic = j-i+2num_w-1; id = i+j+2num_w-3
        r1 = fm*I1[ia] + I2[ia]
        r2 = fm*I1[ib] + I2[ib]
        r3 = fp*I1[ic] + I2[ic]
        r4 = fp*I1[id] + I2[id]

        Kp_real[i,j] =  r1 + r2 - r3 - r4
        Km_real[i,j] =  r1 - r2 + r3 - r4
    end

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
    w_dynam = -(inp.reOmega_c .+ inp.reOmega_c .* cos.((2 .* (inp.n_cheb/2:inp.n_cheb-1) .+ 1) .* π ./ (2 * inp.n_cheb)))    # only positvie chebyshev nodes
    append!(w_dynam, reverse(-w_dynam))   # symmetric chebyshev nodes
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

