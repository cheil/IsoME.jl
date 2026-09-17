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
"""


###########################################
# ------------ Kernel set up ------------ #
###########################################
# Runs once per temperature (entry point `precompute`) and ends in a `PhononKernel`,
# the A(x), B(x) tables the consumer side below evaluates.


"""
    OMEGA_MIN

Lower bound for the first point of the ω-grid. Z(ω) = 1 − Iz(ω)/ω diverges as 1/ω at finite
temperature, so the grid must not start arbitrarily close to zero.
"""
const OMEGA_MIN = 0.1

"""
    omega_grid_start(domega)

First point of the ω-grid: the smallest multiple of `domega` that is at least `OMEGA_MIN`.

Being an *integer multiple* of the grid step is what makes the kernel lattice work out. The
kernel tables are stored on x = x_min + (p−1)·domega with x_min/domega ∈ ℤ, so both folded
arguments of a grid point,

    x_dif = ω' − ω_i,    x_sum = −ω' − ω_i,

land on that same lattice for every row i, and the head lattice of `HeadLattice` is simply
domega·ℤ. With the old hard-coded start of 1.0 this held only for domega ∈ {2, 1, 0.5, 0.4,
0.25, …}; it now holds for every domega. Equals 1.0 at the default domega = 1.0.
"""
omega_grid_start(domega::Real) = ceil(OMEGA_MIN / domega) * domega


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
`domega`, spanning the largest output grid), and the Ω-integration axis `int_axis`.
"""
function setUpOmegaAxis(inp, matval)
    (a2f_omega, a2f) = matval

    # nonzero values of a2F, ensure endpoints are zero. Guarded so an α²F that never rises
    # above the 1e-2 threshold gives a message instead of a MethodError on `nothing` below
    idx_a2f = findfirst(a2f .> 1e-2)                                               # a2F(w>domega) DISCUSS THIS CUTOFF!!
    idx_dw  = findfirst(a2f_omega .> 2*inp.domega)
    isnothing(idx_a2f) && error("α²F stays below 1e-2 over the whole frequency range - check the a2F-file and the selected smearing (ind_smear).\n\n")
    isnothing(idx_dw) && error("The a2F frequency grid does not reach 2·domega = " * string(2*inp.domega) * " meV. Reduce domega or check the a2F-file.\n\n")

    idx_left = min(idx_a2f, idx_dw)
    idx_right = findlast(a2f .> 1e-2)::Int   # a first entry above the threshold implies a last one

    W_left = a2f_omega[idx_left]
    W_right = a2f_omega[idx_right]

    # interpolate
    G = scale(interpolate(a2f[idx_left:idx_right], BSpline(Linear())), a2f_omega[idx_left:idx_right])

    # Frequency and integration grids
    # w_axis must be larger than w_static. Its step is `domega`, the same step the ω- and
    # ω'-grids use: the kernel tables are only ever *sampled* at multiples of domega (the
    # tail is a domega grid, the head lattice below is domega·ℤ), so tabulating them more
    # finely refines a table that is then subsampled at stride domega.
    # All channels (Z, φ, χ) share omega_c, so cDOS and vDOS size w_axis identically.
    w_axis = 0:inp.domega:(inp.omega_c+inp.domega)   # w-axis at which K(ω,ω') is calculated
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
    lo_supp = findfirst(>(0.0), int_axis)
    hi_supp = findlast(<(W_cut), int_axis)

    # an Ω-axis that misses the support of G entirely would only surface as a MethodError
    # in the range below
    isnothing(lo_supp) && error("The Ω-integration axis has no point above 0; the support of α²F is not sampled.\n\n")
    isnothing(hi_supp) && error("The Ω-integration axis has no point below the phonon cutoff ($(W_cut) meV); the support of α²F is not sampled.\n\n")

    idx_supp = lo_supp:hi_supp
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

    ω1 = omega_grid_start(inp.domega)
    # One ω-grid for every channel: Z(ω), Δ(ω)/φ(ω) and χ(ω) are all solved here.
    # Sensitive to the start value, do not choose < 1e-1.
    w_static = ω1:inp.domega:inp.omega_c

    # The kernel tables are only defined on w_axis: the output grid must fit inside it.
    @assert first(w_axis) <= first(w_static) && last(w_static) <= last(w_axis)

    return (w_static, kernel)

end



###########################################
# - Kernel evaluation -- consumer side -- #
###########################################
# Used by wprimeGrid.jl in every Eliashberg iteration: the head is contracted through its
# hat lattice, the tail is summed from the A,B samples. Neither materializes a kernel block.

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



###########################################
# ------------- Head lattice ------------ #
###########################################
# The head grid is a Chebyshev cluster placed on the *integrand's* pole, with spacings down
# to ~1e-5 meV. The kernel it samples there is piecewise linear on the domega lattice, so
# ~10^5 head nodes read the same linear segment: evaluating K at each of them, for every
# output row, would rebuild the same segment over and over.
#
# Contracting over the ω'-nodes first removes that redundancy. With Λ_m the hat functions of
# the head lattice ν_m = (m-1)·domega, the interpolant obeys the identity
#
#     A_h(±ω' − ω_i) = Σ_m Λ_m(ω')·A_h(±ν_m − ω_i),
#
# whose coefficients Λ_m(ω') do not depend on ω_i, and whose values A_h(±ν_m − ω_i) are bare
# table entries (no interpolation) because ±ν_m − ω_i is an exact multiple of domega. Putting
# that into Σ_j w_j g_j 𝒦(ω_i,ω'_j) and exchanging the two finite sums gives
#
#     Σ_j w_j g_j 𝒦(ω_i,ω'_j) = Σ_m Â_{p(m,i)}·G⁰_m + Σ_m B̂_{p(m,i)}·G¹_m,
#     G⁰_m = Σ_j w_j g_j Λ_m(ω'_j),   G¹_m = Σ_j w_j g_j f(ω'_j) Λ_m(ω'_j).
#
# No approximation is involved: the only quadrature is the trapezoidal sum the dense route
# already performs, and the rest is an exchange of finite sums. The Fermi factor is kept in
# its own vector G¹ rather than folded into the hats, so f is never interpolated.
#
# Consequences: no O(n_out × n_head) block is ever formed, the sum runs over
# M ≈ ω'_max/domega instead of n_head, and a moving pole changes only G⁰/G¹ - the kernel is
# untouched, so `maybe_refresh_head!` does no kernel work at all.

"""
    HeadLattice

Head-lattice bookkeeping for the factorised route: the hat lattice ν_m = (m−1)·domega,
m = 1…`M`, covering [0, ω'_max], together with per-head-node data.

  * `p0`, `q0` - the affine table index maps, `p(m,i) = p0 + m − i` (difference channel,
    Toeplitz) and `q(m,i) = q0 − m − i` (sum channel, Hankel). Both are exact integer
    indices into the A/B tables, so the accumulation never interpolates.
  * `cell`, `frac` - the lattice cell m_j and fractional position fr_j of each head node,
    i.e. the same two numbers `_lerpAB` computes, used here to *scatter* into the lattice
    instead of gathering from it.
  * `fp` - f(ω') on the head nodes, evaluated exactly (never interpolated).

Rebuilt whenever the head grid is, which is O(n_head) and involves no kernel data.
"""
struct HeadLattice
    M::Int
    p0::Int
    q0::Int
    cell::Vector{Int}
    frac::Vector{Float64}
    fp::Vector{Float64}
end


"""
    HeadLattice(ker, w_static, wp_head, n_out_max)

Build the lattice for the head grid `wp_head`. `n_out_max` is the largest number of output
rows any channel will ask for (the master-grid length), used to bound-check the table
indices once here so the hot loops can index bare.
"""
function HeadLattice(ker::PhononKernel, w_static, wp_head::Vector{Float64}, n_out_max::Int)
    length(w_static) >= 2 || error("HeadLattice: w_static needs at least two points")
    ω1 = float(w_static[1])
    d = float(w_static[2]) - ω1

    abs(d - ker.dx) <= 1e-9 * max(1.0, d) ||
        error("HeadLattice: ω-grid step $d and kernel table step $(ker.dx) must agree")

    # ν_m = (m−1)·d covering [0, ω'_max]; the +2 leaves room for the upper hat of a node
    # sitting exactly on the last lattice point.
    wp_max = wp_head[end]
    M = floor(Int, wp_max / d) + 2

    # p(m,i) = p0 + m − i and q(m,i) = q0 − m − i follow from the table index of x = ν_1 − ω_1,
    #     t0 = (ν_1 − ω1 − x_min)/dx + 1 = (−ω1 − x_min)/dx + 1,
    # which serves both channels because −ν_1 = ν_1 = 0. That t0 comes out integral IS the
    # lattice condition: it needs ω1/dx ∈ ℤ (guaranteed by `omega_grid_start`) together with
    # x_min/dx ∈ ℤ (guaranteed by `build_kernel`). The check below is therefore the single
    # place where a grid that cannot support the factorisation is caught.
    t0 = (-ω1 - ker.x_min) / ker.dx + 1
    p0 = round(Int, t0)
    abs(t0 - p0) <= 1e-6 || error(
        "HeadLattice: the ω-grid is not commensurate with the kernel table (t0 = $t0). " *
        "The grid start must be an integer multiple of the step - see omega_grid_start.")
    q0 = p0 + 2

    # index ranges the accumulation will touch; checked once so the loops can skip clamping
    pmin, pmax = p0 + 1 - n_out_max, p0 + M - 1
    qmin, qmax = q0 - M - n_out_max, q0 - 2
    (1 <= pmin && pmax <= ker.n_x && 1 <= qmin && qmax <= ker.n_x) || error(
        "HeadLattice: kernel table index out of range (p ∈ [$pmin, $pmax], " *
        "q ∈ [$qmin, $qmax], table has $(ker.n_x) entries)")

    nh = length(wp_head)
    cell = Vector{Int}(undef, nh)
    frac = Vector{Float64}(undef, nh)
    @inbounds for j in 1:nh
        u = wp_head[j] / d
        c = floor(Int, u)
        cell[j] = c + 1
        frac[j] = u - c
    end

    fp = [1 / (exp(ker.β * w) + 1) for w in wp_head]

    return HeadLattice(M, p0, q0, cell, frac, fp)
end


"""
    head_scatter!(G0, G1, lat, wgt, g)

Build the two lattice weight vectors from one integrand density:

    G⁰_m = Σ_j w_j g_j Λ_m(ω'_j),    G¹_m = Σ_j w_j g_j f(ω'_j) Λ_m(ω'_j).

The adjoint of linear interpolation: each head node deposits into exactly two lattice cells,
with the weights (1−fr, fr) that `_lerpAB` would have used to read out of them. O(n_head),
independent of the number of output rows.
"""
function head_scatter!(G0::Vector{Float64}, G1::Vector{Float64}, lat::HeadLattice,
                       wgt::AbstractVector{Float64}, g::AbstractVector{Float64})
    fill!(G0, 0.0)
    fill!(G1, 0.0)
    @inbounds for j in eachindex(wgt)
        h = wgt[j] * g[j]
        m = lat.cell[j]
        fr = lat.frac[j]
        w_lo = h * (1 - fr)
        w_hi = h * fr
        f = lat.fp[j]
        G0[m]   += w_lo
        G0[m+1] += w_hi
        G1[m]   += w_lo * f
        G1[m+1] += w_hi * f
    end
    return
end


"""
    head_accumulate!(Iz, Iphi, ker, lat, G0z, G1z, G0p, G1p, nout)

Head contribution to the K⁻ (Iz) and K⁺ (Iphi) channels, *overwriting* the first `nout`
entries. With
𝒦_dif = Â_p + f·B̂_p and 𝒦_sum = Â_q + (1−f)·B̂_q, so that K⁻ = −𝒦_dif − 𝒦_sum and
K⁺ = 𝒦_dif − 𝒦_sum, the (1−f) of the sum channel is what splits B̂_q across both weight
vectors, with opposite relative signs in the two channels.
"""
function head_accumulate!(Iz::Vector{ComplexF64}, Iphi::Vector{ComplexF64},
                          ker::PhononKernel, lat::HeadLattice,
                          G0z::Vector{Float64}, G1z::Vector{Float64},
                          G0p::Vector{Float64}, G1p::Vector{Float64}, nout::Int)
    A = ker.A
    B = ker.B
    M = lat.M
    @inbounds for i in 1:nout
        accz = zero(ComplexF64)
        accp = zero(ComplexF64)
        p = lat.p0 + 1 - i          # p(1,i), then +1 per lattice cell
        q = lat.q0 - 1 - i          # q(1,i), then -1 per lattice cell
        for m in 1:M
            Ap = A[p]; Bp = B[p]
            Aq = A[q]; Bq = B[q]
            accz += (Ap + Aq + Bq) * G0z[m] + (Bp - Bq) * G1z[m]
            accp += (Ap - Aq - Bq) * G0p[m] + (Bp + Bq) * G1p[m]
            p += 1
            q -= 1
        end
        Iz[i] = -accz
        Iphi[i] = accp
    end
    return
end


"""
    head_accumulate3!(Iz, Iphi, Ichi, ker, lat, G0z, G1z, G0p, G1p, G0c, G1c, nout)

Three-channel head contribution, *overwriting* the first `nout` entries of all three.

χ uses the same K⁺ combination as φ, so it rides along on the A/B lookups that Iz and Iphi
already perform: one pass over the kernel tables instead of the two it used to take when χ
lived on its own (wider) grid.
"""
function head_accumulate3!(Iz::Vector{ComplexF64}, Iphi::Vector{ComplexF64}, Ichi::Vector{ComplexF64},
                           ker::PhononKernel, lat::HeadLattice,
                           G0z::Vector{Float64}, G1z::Vector{Float64},
                           G0p::Vector{Float64}, G1p::Vector{Float64},
                           G0c::Vector{Float64}, G1c::Vector{Float64}, nout::Int)
    A = ker.A
    B = ker.B
    M = lat.M
    @inbounds for i in 1:nout
        accz = zero(ComplexF64)
        accp = zero(ComplexF64)
        accc = zero(ComplexF64)
        p = lat.p0 + 1 - i          # p(1,i), then +1 per lattice cell
        q = lat.q0 - 1 - i          # q(1,i), then -1 per lattice cell
        for m in 1:M
            Ap = A[p]; Bp = B[p]
            Aq = A[q]; Bq = B[q]
            Kp0 = Ap - Aq - Bq      # K⁺ weight of G0, shared by the φ and χ channels
            Kp1 = Bp + Bq           # K⁺ weight of G1
            accz += (Ap + Aq + Bq) * G0z[m] + (Bp - Bq) * G1z[m]
            accp += Kp0 * G0p[m] + Kp1 * G1p[m]
            accc += Kp0 * G0c[m] + Kp1 * G1c[m]
            p += 1
            q -= 1
        end
        Iz[i] = -accz
        Iphi[i] = accp
        Ichi[i] = accc
    end
    return
end


"""
    tail_accumulate3!(Iz, Iphi, Ichi, A_dif, B_dif, A_sum, B_sum, noff, nout, wgt, fp, fm,
                      g_z, g_phi, g_chi)

Three-channel tail accumulation. Same indexing as `tail_accumulate!`; χ reuses the K⁺
already formed for φ, so the A/B lookups are done once for all three channels.
"""
function tail_accumulate3!(Iz::Vector{ComplexF64}, Iphi::Vector{ComplexF64}, Ichi::Vector{ComplexF64},
                           A_dif::Vector{ComplexF64}, B_dif::Vector{ComplexF64},
                           A_sum::Vector{ComplexF64}, B_sum::Vector{ComplexF64},
                           noff::Int, nout::Int,
                           wgt::AbstractVector{Float64},
                           fp::AbstractVector{Float64}, fm::AbstractVector{Float64},
                           g_z::AbstractVector{Float64}, g_phi::AbstractVector{Float64},
                           g_chi::AbstractVector{Float64})
    nt = length(wgt)
    @inbounds for j in 1:nt
        fpj = fp[j]; fmj = fm[j]
        hz = wgt[j] * g_z[j]
        hp = wgt[j] * g_phi[j]
        hc = wgt[j] * g_chi[j]
        base_dif = j + noff      # p = base_dif - i
        base_sum = j             # q = base_sum + i - 1

        for i in 1:nout
            p = base_dif - i
            q = base_sum + i - 1
            K_dif = A_dif[p] + fpj * B_dif[p]     # 𝒦(ω, ω')
            K_sum = A_sum[q] + fmj * B_sum[q]     # 𝒦(ω,-ω')
            Kplus = K_dif - K_sum                 # K⁺, shared by the φ and χ channels

            Iz[i]   += (-K_dif - K_sum) * hz      # K⁻
            Iphi[i] += Kplus * hp
            Ichi[i] += Kplus * hc
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

# The in-place bar redraws itself with "\r", which only blanks the line on a real terminal.
# With stdout redirected (batch job, CI log, `julia run.jl > out.txt`) every frame would
# instead be appended, so the bar is drawn only on a TTY. The completed bar written to the
# log file by `finish_kernel_progress` is unaffected.
isProgressTTY() = stdout isa Base.TTY

function redraw_kernel_progress(progress_state)
    spinner = raw"-\|/"[mod(progress_state.spinner[], 4) + 1]
    progress_state.spinner[] += 1
    isProgressTTY() || return
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
    # off a TTY nothing was drawn, so print the finished bar once instead of a bare newline
    isProgressTTY() || print(progress_state.prefix * "|" * "="^progress_state.target * "|")
    print("\n")
    if !isnothing(log_file)
        print(log_file, progress_state.prefix * "|" *  "="^progress_state.target * "|\n")
    end
end
