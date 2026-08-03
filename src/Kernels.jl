"""
    Real-axis phonon kernels: build and evaluation.

    The kernel is never stored as a matrix. K(ω,ω') is built from two *one-dimensional*
    complex functions of the difference variable x, times the Fermi factor of ω':

        𝒦(ω,ω') = A(x) + f(ω')·B(x),        x = ω' - ω

        K⁺(ω,ω') = +𝒦(ω,ω') - 𝒦(ω,-ω')
        K⁻(ω,ω') = -𝒦(ω,ω') - 𝒦(ω,-ω')

    with A, B assembled from the principal-value integrals over α²F,

        I1(x) = P∫dΩ G(Ω)/(Ω-x),    I2(x) = P∫dΩ G(Ω)·n(Ω)/(Ω-x)

        A(x) = I1(-x) + I2(-x) - I2(x) + iπ[ G(x)·n(x) + G(-x)·(1 + n(-x)) ]
        B(x) = -I1(x) - I1(-x)         + iπ[ G(x) - G(-x) ]

    Only the O(N) tables A, B are kept; the ω'-integral evaluates A + f·B on the fly.

    SET UP SIDE (once per temperature, entry point `precompute`)

        precompute -> setUpOmegaAxis          α²F -> G, w_axis, int_axis
                   -> kernels                 -> compute_kernel_integrals
                                                 -> precompute_integrals
                                                    -> kernel_integral_helper   I1, I2
                                              -> build_kernel -> assemble_AB    A, B
                                                 => PhononKernel

    CONSUMER SIDE (every Eliashberg iteration, driven by wprimeGrid.jl)

        WPrimeWorkspace  -> fermi_pm, precompute_tail_AB       tail samples, once per grid
                         -> fill_kp!, fill_km!                 materialized head blocks
        ω'-integral      -> tail_accumulate!, tail_accumulate_chi!
"""


###########################################
# ------------ Kernel set up ------------ #
###########################################
# Runs once per temperature (entry point `precompute`) and ends in a `PhononKernel`,
# the A(x), B(x) tables the consumer side below evaluates.


"""
    PhononKernel

The kernel as the two tables A(x), B(x) on a uniform x-grid (`n_x` points from `x_min`, step
`dx`) plus the inverse temperature β it was built at. O(N) storage; no 2-D matrix.
"""
struct PhononKernel
    x_min::Float64
    dx::Float64
    n_x::Int
    A::Vector{ComplexF64}
    B::Vector{ComplexF64}
    β::Float64
end


"""
    assemble_AB(xgrid, I1, I2, G, W_left, W_right, β)

A and B from the principal-value integrals I1, I2 evaluated on `xgrid` (which must be
symmetric, so that -x is the mirrored index).
"""
function assemble_AB(xgrid::Vector{Float64}, I1::Vector{Float64}, I2::Vector{Float64},
                     G, W_left::Float64, W_right::Float64, β::Float64)
    L = length(xgrid)
    nb = y -> 1 / (exp(β * y) - 1)      # Bose function, only called for y > W_left > 0

    A = Vector{ComplexF64}(undef, L)
    B = Vector{ComplexF64}(undef, L)

    @inbounds for p in 1:L
        x = xgrid[p]
        q = L + 1 - p                       # index of -x (grid is symmetric)

        i1x = I1[p]; i1m = I1[q]
        i2x = I2[p]; i2m = I2[q]

        # at most one of the two is nonzero (G lives on [W_left, W_right], W_left > 0)
        inx = W_left < x < W_right
        inm = W_left < -x < W_right
        Gx = inx ? G(x) : 0.0
        Gm = inm ? G(-x) : 0.0

        # on-shell (δ-function) part; the Bose argument is positive wherever G ≠ 0
        imA = inx ? Gx * nb(x) : (inm ? Gm * (1 + nb(-x)) : 0.0)

        A[p] = (i1m + i2m - i2x) + im * π * imA
        B[p] = -(i1x + i1m) + im * π * (Gx - Gm)
    end

    return A, B
end


"""
    build_kernel(β, num_w, w_max, I1, I2, G, W_left, W_right)

Assemble the A, B tables on the uniform difference grid `precompute_integrals` produced
(step dΩ, reaching ±2·w_max). O(N) storage.
"""
function build_kernel(β::Float64, num_w::Int, w_max::Float64,
                             I1::Vector{Float64}, I2::Vector{Float64},
                             G, W_left::Float64, W_right::Float64)

    L = 4 * num_w - 3
    length(I1) == L || error("I1 has length $(length(I1)), expected $L")

    dx = w_max / (num_w - 1)
    X_max = 2 * (num_w - 1) * dx
    xgrid = collect(range(-X_max, X_max, length=L))

    A, B = assemble_AB(xgrid, I1, I2, G, W_left, W_right, β)

    return PhononKernel(-X_max, dx, L, A, B, β)
end


"""
    make_integration_axis(W_width, pts_cheb, pts_lin, epsilon)

Set up the Ω-integration axis: a Chebyshev grid inside [-epsilon, epsilon], where the 1/Ω
pole sits, and a linear grid out to ±W_width on either side.
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
    setUpOmegaAxis(inp, matval)

Everything the kernel build needs from α²F: the interpolant `G` and its support
[`W_left`, `W_right`], the ω-grid `w_axis` on which the kernel tables are laid out (step
`dKernel`, spanning the largest output grid), and the Ω-integration axis `int_axis`.
"""
function setUpOmegaAxis(inp, matval)
    (a2f_omega, a2f) = matval

    # nonzero values of a2F, ensure endpoints are zero
    idx_left = min(findfirst(a2f .> 1e-2), findfirst(a2f_omega .> 2*inp.domega))   # a2F(w>domega) DISCUSS THIS CUTOFF!!
    idx_right = findlast(a2f .> 1e-2)

    W_left = a2f_omega[idx_left]
    W_right = a2f_omega[idx_right]

    # interpolate
    G = scale(interpolate(a2f[idx_left:idx_right], BSpline(Linear())), a2f_omega[idx_left:idx_right])

    # Frequency and integration grids
    # w_axis must be larger than max(w_static, w_static_chi)
    if inp.cDOS_flag == 1
        w_axis = 0:inp.dKernel:(inp.reOmega_c+inp.dKernel)   # w-axis at which K(ω,ω') is calculated
    else
        # ω linear grid to store K(ω,ω'), has to span max(reOmega_c, reOmega_c_shift).
        w_axis = 0:inp.dKernel:(max(inp.reOmega_c, inp.reOmega_c_shift)+inp.dKernel)
    end
    # Ω-integration axis: 300 Chebyshev points inside ±3 meV (the 1/Ω pole), 300 linear
    # points out to ±(W_right-W_left), the width over which G is nonzero.
    int_axis = make_integration_axis(W_right-W_left, 300, 300, 3)

    return G, W_left, W_right, w_axis, int_axis
end


"""
    kernel_integral_helper(x_grid, int_axis, W_cut, W_left, G, bose; progress...)

The two principal-value integrals on `x_grid`, by trapezoidal Ω-integration over `int_axis`:

    I1(x) = P∫dΩ G(Ω)/(Ω-x),    I2(x) = P∫dΩ G(Ω)·n(Ω)/(Ω-x)

Two cases, distinguished by whether the pole at Ω = x falls inside the support of G:

  * `0 < x < W_cut` — the pole sits inside the integration range. The integrand is written
    around it, G(Ω+x)/Ω, and has to be evaluated point by point. Only a handful of x are of
    this kind (as many as the α²F bandwidth holds, ~200 of ~100000).
  * otherwise — the pole is outside, the integrand is G(Ω)/(Ω-x) and factorizes into a part
    that depends on Ω alone and the factor 1/(Ω-x). The Ω-only parts (`wG`, `wGn`) are
    tabulated once and every such x is then a dot product with 1/(Ω-x).

Since G vanishes outside (0, W_cut) only that contiguous stretch of `int_axis` contributes;
dropping the exact zeros elsewhere leaves the trapezoidal sum unchanged.

Returns (I1, I2), each of length `length(x_grid)`.
"""
function kernel_integral_helper(
    x_grid::AbstractVector{Float64},
    int_axis::Vector{Float64},
    W_cut::Float64,
    W_left::Float64,
    G::AbstractInterpolation,
    bose;
    progress_counter=nothing,
    progress_interval=nothing,
    progress_state=nothing,
)

    n_x = length(x_grid)
    n_Ω = length(int_axis)

    # trapezoid weights of the Ω-axis: Σ_j wgt_Ω[j]·f(Ω_j) == trapz(int_axis, f)
    wgt_Ω = similar(int_axis)
    wgt_Ω[1]   = (int_axis[2] - int_axis[1]) / 2
    wgt_Ω[n_Ω] = (int_axis[n_Ω] - int_axis[n_Ω-1]) / 2
    @inbounds for j in 2:n_Ω-1
        wgt_Ω[j] = (int_axis[j+1] - int_axis[j-1]) / 2
    end

    # support of G on the Ω-axis, and the Ω-only parts of the integrands there
    idx_supp = findfirst(>(0.0), int_axis):findlast(<(W_cut), int_axis)
    Ω_supp   = collect(int_axis[idx_supp])
    wG       = [wgt_Ω[j] * G(int_axis[j] + W_left) for j in idx_supp]
    wGn      = [wG[k] * bose(Ω_supp[k] + W_left) for k in eachindex(Ω_supp)]

    integral_1 = zeros(n_x)
    integral_2 = zeros(n_x)
    @inbounds for i in 1:n_x
        x = x_grid[i]
        acc_1 = 0.0
        acc_2 = 0.0

        if 0 < x < W_cut
            # pole inside the support: evaluate G and n point by point
            for j in 1:n_Ω
                Ω = int_axis[j]
                if 0 < Ω + x < W_cut
                    integrand_1 = G(Ω + x + W_left) / Ω
                    acc_1 += wgt_Ω[j] * integrand_1
                    acc_2 += wgt_Ω[j] * integrand_1 * bose(Ω + x + W_left)
                end
            end
        else
            # pole outside: only 1/(Ω-x) depends on x, the rest is tabulated
            @simd for k in eachindex(Ω_supp)
                pole_factor = 1 / (Ω_supp[k] - x)
                acc_1 += wG[k]  * pole_factor
                acc_2 += wGn[k] * pole_factor
            end
        end

        integral_1[i] = acc_1
        integral_2[i] = acc_2

        if !isnothing(progress_counter) && !isnothing(progress_interval)
            progress_counter[] += 1
            if mod(progress_counter[], progress_interval) == 0
                print_kernel_progress(progress_state)
            end
        end
    end

    return integral_1, integral_2
end


"""
    precompute_integrals(num_w, w_max, W_left, W_right, int_axis, G, bose; progress_state)

I1, I2 on the unified difference grid x_m = m·dw, m = -2(num_w-1) … 2(num_w-1) — every
value the ±ω±ω' combinations of the kernel can reach, so a single 1-D table serves all four
folded arguments. Length 4·num_w-3.
"""
function precompute_integrals(num_w, w_max, W_left, W_right, int_axis, G, bose; progress_state=nothing)
    W_cut = W_right - W_left
    dw = w_max / (num_w - 1)
    progress_counter = isnothing(progress_state) ? nothing : Ref(0)
    progress_interval = isnothing(progress_counter) ? nothing : kernel_progress_interval(3 * (2 * num_w - 1), progress_state)

    # single unified x-grid: position p ↔ frequency index m = p - (2*num_w - 1),
    # covering every m ∈ -2(num_w-1) … 2(num_w-1) that the ±w±w' combinations reach.
    x_idx  = -2*(num_w-1):2*(num_w-1)
    x_grid = x_idx .* dw .- W_left
    return kernel_integral_helper(
        x_grid,
        int_axis,
        W_cut,
        W_left,
        G,
        bose;
        progress_counter=progress_counter,
        progress_interval=progress_interval,
        progress_state=progress_state,
    )
end

"""
    compute_kernel_integrals(num_w, w_max, W_left, W_right, int_axis, G, bose, log_file, line_width)

`precompute_integrals` plus the progress bar. I1(x)=P∫G(Ω)/(Ω−x) and I2(x)=P∫G(Ω)n(Ω)/(Ω−x)
are the only kernel data the kernel needs (𝒦 = A(x) + f(ω')·B(x)); this is by far the
most expensive step of the setup.
"""
function compute_kernel_integrals(num_w, w_max, W_left, W_right, int_axis, G, bose, log_file=nothing, line_width=80)

    progress_state = isnothing(log_file) ? nothing : start_kernel_progress(" ", line_width)

    I1, I2 = precompute_integrals(num_w, w_max, W_left, W_right, int_axis, G, bose; progress_state=progress_state)

    isnothing(log_file) || finish_kernel_progress(log_file, progress_state)

    return I1, I2
end


"""
    kernels(β, w_axis, W_left, W_right, int_axis, G, log_file, line_width)

Build the kernel at a given temperature: the 1-D tables A(x), B(x) with
    𝒦(ω,ω') = A(x) + f(ω')·B(x),   x = ω'−ω
    Kp(ω,ω') = +𝒦(ω,ω') − 𝒦(ω,-ω')
    Km(ω,ω') = -𝒦(ω,ω') − 𝒦(ω,-ω')
with K(ω,ω') as given in the paper, supplemental eq. (47). Returns the `PhononKernel`; the
dense O(N²) Kp/Km matrices are never formed. β enters twice: through the Bose factor of the
Ω-integral, and through the Fermi factor assembled into A, B.
"""
function kernels(β, w_axis, W_left, W_right, int_axis, G, log_file=nothing, line_width=80)

    # Bose-Einstein occupation of the phonons
    bose = x -> x == 0 ? 0.0 : abs(1 / (exp(β * x) - 1))

    I1, I2 = compute_kernel_integrals(length(w_axis), Float64(w_axis[end]), W_left, W_right,
                                      int_axis, G, bose, log_file, line_width)

    return build_kernel(β, length(w_axis), Float64(w_axis[end]), I1, I2, G, W_left, W_right)
end

function precompute(β, inp, matval, console, log_file)

    printTextCentered("Precomputing Kernels", console.cDOS.partingLine, file=log_file, bold=true)
    print_kernel_log(log_file, "\n")

    ### Set up axis and parameters ###
    G, W_left, W_right, w_axis, int_axis = setUpOmegaAxis(inp, matval)

    kernel = kernels(β, w_axis, W_left, W_right, int_axis, G, log_file, length(console.cDOS.partingLine))

    w_static = 1:inp.domega:inp.reOmega_c        # grid of Z(w), Delta(w), sensitive to start value, do not chose < 1e-1
    w_static_chi = 1:inp.domega:inp.reOmega_c_shift   # grid of χ(ω); same step as w_static so all channels are k=1

    # The kernel tables are only defined on w_axis: both output grids must fit inside it.
    # w_static_chi exists in every vDOS mode (χ is solved there), not just vDOS+W, and
    # setUpOmegaAxis sizes w_axis accordingly.
    @assert first(w_axis) <= first(w_static) && last(w_static) <= last(w_axis)
    if inp.cDOS_flag == 0
        @assert first(w_axis) <= first(w_static_chi) && last(w_static_chi) <= last(w_axis)
    end

    return (w_static, w_static_chi, kernel)

end



###########################################
# - Kernel evaluation -- consumer side -- #
###########################################
# Used by wprimeGrid.jl in every Eliashberg iteration: the head blocks are
# materialized once per ω'-grid, the tail is summed from the A,B samples.

"""
    fermi_pm(β, wp)

Precompute the Fermi factors f(ω') and f(-ω') on the ω'-grid `wp`. O(N).
"""
function fermi_pm(β::Float64, wp::AbstractVector{<:Real})
    fp = [1 / (exp(β * w) + 1) for w in wp]
    fm = 1 .- fp
    return fp, fm
end


"""
    _lerpAB(A, B, t, n_x)

Linear interpolation of A and B at the fractional (1-based) grid index t, clamped to the
grid. Sharing t between the two tables halves the index arithmetic.
"""
@inline function _lerpAB(A::Vector{ComplexF64}, B::Vector{ComplexF64}, t::Float64, n_x::Int)
    k = unsafe_trunc(Int, t)
    k = ifelse(k < 1, 1, ifelse(k > n_x - 1, n_x - 1, k))
    fr = t - k
    @inbounds begin
        a = A[k]; b = B[k]
        return (a + fr * (A[k+1] - a), b + fr * (B[k+1] - b))
    end
end


"""
    kernel_eval(ker, w, wp)

K⁺(ω,ω') and K⁻(ω,ω') straight from the A, B tables (single-point convenience/debug).
"""
@inline function kernel_eval(ker::PhononKernel, w::Float64, wp::Float64)
    invdx = 1 / ker.dx
    A_dif, B_dif = _lerpAB(ker.A, ker.B, (wp - w - ker.x_min) * invdx + 1, ker.n_x)   # x = ω'-ω
    A_sum, B_sum = _lerpAB(ker.A, ker.B, (-wp - w - ker.x_min) * invdx + 1, ker.n_x)  # x = -ω-ω'

    fp = 1 / (exp(ker.β * wp) + 1)        # f(ω')
    fm = 1 - fp                          # f(-ω')

    K_dif = A_dif + fp * B_dif           # 𝒦(ω, ω')
    K_sum = A_sum + fm * B_sum           # 𝒦(ω,-ω')

    return (K_dif - K_sum, -K_dif - K_sum)   # (K⁺, K⁻)
end


"""
    precompute_tail_AB(ker, w_static, wp_tail)

Sample A, B onto the two folded arguments of the *uniform* tail grid:

    difference channel   x_dif = ω'_j - ω_i  =  c_dif + (j-i)·dw   → builds 𝒦(ω, ω')
    sum channel          x_sum = -ω_i - ω'_j =  c_sum - (i+j-2)·dw → builds 𝒦(ω,-ω')

Because w_static and wp_tail share the step dw, each argument depends on a single offset
(j-i for the difference / Toeplitz, i+j for the sum / Hankel), so only O(ns + n_tail)
samples are stored instead of the ns×n_tail matrix. Returns (A_dif, B_dif, A_sum, B_sum),
each length ns + n_tail - 1: A_dif/B_dif indexed by p = (j-i)+ns, A_sum/B_sum by q = i+j-1.
The A,B interpolation is done once here; the ω'-integral then only looks them up.
"""
function precompute_tail_AB(ker::PhononKernel, w_static, wp_tail::AbstractVector{Float64})
    ns = length(w_static)
    nt = length(wp_tail)
    invdx = 1 / ker.dx
    xmin = ker.x_min

    ws1 = float(w_static[1])
    wt1 = float(wp_tail[1])
    dw = float(w_static[2]) - ws1
    abs((wp_tail[2] - wt1) - dw) < 1e-9 * max(1.0, dw) ||
        error("precompute_tail_AB: w_static and wp_tail must share the same step")

    L = ns + nt - 1
    A_dif = Vector{ComplexF64}(undef, L); B_dif = Vector{ComplexF64}(undef, L)
    A_sum = Vector{ComplexF64}(undef, L); B_sum = Vector{ComplexF64}(undef, L)

    c_dif = wt1 - ws1          # x_dif(p) = c_dif + (p - ns)·dw
    c_sum = -(wt1 + ws1)       # x_sum(q) = c_sum - (q - 1)·dw
    @inbounds for p in 1:L
        x = c_dif + (p - ns) * dw
        A_dif[p], B_dif[p] = _lerpAB(ker.A, ker.B, (x - xmin) * invdx + 1, ker.n_x)
    end
    @inbounds for q in 1:L
        x = c_sum - (q - 1) * dw
        A_sum[q], B_sum[q] = _lerpAB(ker.A, ker.B, (x - xmin) * invdx + 1, ker.n_x)
    end
    return A_dif, B_dif, A_sum, B_sum
end


"""
    fill_kp!(Kp, w, wp, ker)   /   fill_km!(Km, w, wp, ker)

Fill only the K⁺ (resp. K⁻) block on the grid (w, wp) from the A,B tables. K⁺ is built on the
master grid (χ for vDOS), K⁻ only on w_static.
"""
function fill_kp!(Kp::Matrix{ComplexF64}, w, wp::Vector{Float64}, ker::PhononKernel)
    size(Kp) == (length(w), length(wp)) || throw(DimensionMismatch("Kp block does not match the grids"))
    @inbounds for i in eachindex(wp)
        b = wp[i]
        for j in eachindex(w)
            kp, _ = kernel_eval(ker, float(w[j]), b)
            Kp[j, i] = kp
        end
    end
    return Kp
end

function fill_km!(Km::Matrix{ComplexF64}, w, wp::Vector{Float64}, ker::PhononKernel)
    size(Km) == (length(w), length(wp)) || throw(DimensionMismatch("Km block does not match the grids"))
    @inbounds for i in eachindex(wp)
        b = wp[i]
        for j in eachindex(w)
            _, km = kernel_eval(ker, float(w[j]), b)
            Km[j, i] = km
        end
    end
    return Km
end


"""
    tail_accumulate_chi!(Ichi, A_dif, B_dif, A_sum, B_sum, noff, nout, wgt, fp, fm, g_chi)

K⁺-only tail accumulation for the χ channel (Ichi += ∫K⁺ g_chi over the tail). `noff` is the
difference-index offset of the shared tail arrays (their build length = master-grid length),
`nout` the number of output ω-points (here nout = noff = nchi).
"""
function tail_accumulate_chi!(Ichi::Vector{ComplexF64},
                                  A_dif::Vector{ComplexF64}, B_dif::Vector{ComplexF64},
                                  A_sum::Vector{ComplexF64}, B_sum::Vector{ComplexF64},
                                  noff::Int, nout::Int,
                                  wgt::AbstractVector{Float64},
                                  fp::AbstractVector{Float64}, fm::AbstractVector{Float64},
                                  g_chi::AbstractVector{Float64})
    nt = length(wgt)
    @inbounds for j in 1:nt
        fpj = fp[j]; fmj = fm[j]
        hc = wgt[j] * g_chi[j]
        base_dif = j + noff
        base_sum = j

        for i in 1:nout
            p = base_dif - i
            q = base_sum + i - 1
            K_dif = A_dif[p] + fpj * B_dif[p]
            K_sum = A_sum[q] + fmj * B_sum[q]
            Ichi[i] += (K_dif - K_sum) * hc      # K⁺
        end
    end
    return
end


"""
    tail_accumulate!(Iz, Iphi, A_dif, B_dif, A_sum, B_sum, ns, wgt, fp, fm, g_z, g_phi)

Accumulate the tail ω'-segment into Iz/Iphi from the precomputed samples — one loop over ω'
building the K_tail column for all ω by pure lookups (𝒦 = A + f·B), no interpolation.
`A_dif`/`B_dif` are indexed by p = (j-i)+noff, `A_sum`/`B_sum` by q = i+j-1. `noff` is the
difference offset of the shared tail arrays (their build/master length); `nout` is the number
of ω-output points (ns for the Km/Kp channels, which loop only the first nout ≤ noff rows).
"""
function tail_accumulate!(Iz::Vector{ComplexF64}, Iphi::Vector{ComplexF64},
                              A_dif::Vector{ComplexF64}, B_dif::Vector{ComplexF64},
                              A_sum::Vector{ComplexF64}, B_sum::Vector{ComplexF64},
                              noff::Int, nout::Int,
                              wgt::AbstractVector{Float64},
                              fp::AbstractVector{Float64}, fm::AbstractVector{Float64},
                              g_z::AbstractVector{Float64}, g_phi::AbstractVector{Float64})
    nt = length(wgt)
    @inbounds for j in 1:nt
        fpj = fp[j]; fmj = fm[j]
        hz = wgt[j] * g_z[j]
        hp = wgt[j] * g_phi[j]
        base_dif = j + noff      # p = base_dif - i
        base_sum = j             # q = base_sum + i - 1

        for i in 1:nout
            p = base_dif - i
            q = base_sum + i - 1
            K_dif = A_dif[p] + fpj * B_dif[p]     # 𝒦(ω, ω')
            K_sum = A_sum[q] + fmj * B_sum[q]     # 𝒦(ω,-ω')

            Iz[i]   += (-K_dif - K_sum) * hz      # K⁻
            Iphi[i] += (K_dif - K_sum) * hp       # K⁺
        end
    end
    return
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

function start_kernel_progress(label, line_width)
    prefix = label * " K(w,w'): "
    progress_state = (
        printed=Ref(0),
        target=max(1, line_width - length(prefix) - 3),
        spinner=Ref(0),
        prefix=prefix,
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

function print_kernel_progress(progress_state)
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
