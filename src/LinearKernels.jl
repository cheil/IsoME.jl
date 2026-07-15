"""
    Linear (difference-grid) representation of the real-axis kernels.

    K(ω,ω') is not a function of ω'-ω alone, but it is built from two *one-dimensional*
    complex functions of the difference variable x, times the Fermi factor of ω':

        𝒦(ω,ω') = A(x) + f(ω')·B(x),        x = ω' - ω

        A(x) = I1(-x) + I2(-x) - I2(x) + iπ[ G(x)·n(x) + G(-x)·(1 + n(-x)) ]
        B(x) = -I1(x) - I1(-x)          + iπ[ G(x) - G(-x) ]

    with the principal-value integrals `precompute_integrals` already returns,

        I1(x) = P∫dΩ G(Ω)/(Ω-x),    I2(x) = P∫dΩ G(Ω)·n(Ω)/(Ω-x),

    tabulated on the uniform x-grid x_p = x_min + (p-1)·dx, p = 1 … 4·num_w-3.

    The kernels follow by folding ω' → -ω':

        K⁺(ω,ω') = +𝒦(ω,ω') - 𝒦(ω,-ω') = +[A(ω'-ω)+f(ω')B(ω'-ω)] - [A(-ω-ω')+f(-ω')B(-ω-ω')]
        K⁻(ω,ω') = -𝒦(ω,ω') - 𝒦(ω,-ω') = -[A(ω'-ω)+f(ω')B(ω'-ω)] - [A(-ω-ω')+f(-ω')B(-ω-ω')]

    Only the length-O(N) tables A, B (and, per ω'-grid, the Fermi factors f(±ω')) are stored;
    the O(N²) kernel matrix is never formed. The ω'-integral evaluates A + f·B on the fly.
"""


"""
    LinearKernel

The A(x), B(x) tables on a uniform x-grid (x_min, dx, n points) and the inverse temperature β.
Length O(N); no 2-D matrix.
"""
struct LinearKernel
    x_min::Float64
    dx::Float64
    n::Int
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
    build_linear_kernel(β, num_w, w_max, I1, I2, G, W_left, W_right)

Assemble the A, B tables on the uniform difference grid `precompute_integrals` produced
(step dΩ, reaching ±2·w_max). O(N) storage.
"""
function build_linear_kernel(β::Float64, num_w::Int, w_max::Float64,
                             I1::Vector{Float64}, I2::Vector{Float64},
                             G, W_left::Float64, W_right::Float64)

    L = 4 * num_w - 3
    length(I1) == L || error("I1 has length $(length(I1)), expected $L")

    dx = w_max / (num_w - 1)
    X_max = 2 * (num_w - 1) * dx
    xgrid = collect(range(-X_max, X_max, length=L))

    A, B = assemble_AB(xgrid, I1, I2, G, W_left, W_right, β)

    return LinearKernel(-X_max, dx, L, A, B, β)
end


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
    _lerpAB(A, B, t, n)

Linear interpolation of A and B at the fractional (1-based) grid index t, clamped to the
grid. Sharing t between the two tables halves the index arithmetic.
"""
@inline function _lerpAB(A::Vector{ComplexF64}, B::Vector{ComplexF64}, t::Float64, n::Int)
    k = unsafe_trunc(Int, t)
    k = ifelse(k < 1, 1, ifelse(k > n - 1, n - 1, k))
    fr = t - k
    @inbounds begin
        a = A[k]; b = B[k]
        return (a + fr * (A[k+1] - a), b + fr * (B[k+1] - b))
    end
end


"""
    lin_kernel_eval(lk, w, wp)

K⁺(ω,ω') and K⁻(ω,ω') straight from the A, B tables (single-point convenience/debug).
"""
@inline function lin_kernel_eval(lk::LinearKernel, w::Float64, wp::Float64)
    invdx = 1 / lk.dx
    A_dif, B_dif = _lerpAB(lk.A, lk.B, (wp - w - lk.x_min) * invdx + 1, lk.n)   # x = ω'-ω
    A_sum, B_sum = _lerpAB(lk.A, lk.B, (-wp - w - lk.x_min) * invdx + 1, lk.n)  # x = -ω-ω'

    fp = 1 / (exp(lk.β * wp) + 1)        # f(ω')
    fm = 1 - fp                          # f(-ω')

    K_dif = A_dif + fp * B_dif           # 𝒦(ω, ω')
    K_sum = A_sum + fm * B_sum           # 𝒦(ω,-ω')

    return (K_dif - K_sum, -K_dif - K_sum)   # (K⁺, K⁻)
end


"""
    precompute_tail_AB(lk, w_static, wp_tail)

Sample A, B onto the two folded arguments of the *uniform* tail grid:

    difference channel   x_dif = ω'_j - ω_i  =  c_dif + (j-i)·dw   → builds 𝒦(ω, ω')
    sum channel          x_sum = -ω_i - ω'_j =  c_sum - (i+j-2)·dw → builds 𝒦(ω,-ω')

Because w_static and wp_tail share the step dw, each argument depends on a single offset
(j-i for the difference / Toeplitz, i+j for the sum / Hankel), so only O(ns + n_tail)
samples are stored instead of the ns×n_tail matrix. Returns (A_dif, B_dif, A_sum, B_sum),
each length ns + n_tail - 1: A_dif/B_dif indexed by p = (j-i)+ns, A_sum/B_sum by q = i+j-1.
The A,B interpolation is done once here; the ω'-integral then only looks them up.
"""
function precompute_tail_AB(lk::LinearKernel, w_static, wp_tail::AbstractVector{Float64})
    ns = length(w_static)
    nt = length(wp_tail)
    invdx = 1 / lk.dx
    xmin = lk.x_min

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
        A_dif[p], B_dif[p] = _lerpAB(lk.A, lk.B, (x - xmin) * invdx + 1, lk.n)
    end
    @inbounds for q in 1:L
        x = c_sum - (q - 1) * dw
        A_sum[q], B_sum[q] = _lerpAB(lk.A, lk.B, (x - xmin) * invdx + 1, lk.n)
    end
    return A_dif, B_dif, A_sum, B_sum
end


"""
    fill_kp_lin!(Kp, w, wp, lk)   /   fill_km_lin!(Km, w, wp, lk)

Fill only the K⁺ (resp. K⁻) block on the grid (w, wp) from the A,B tables. K⁺ is built on the
master grid (χ for vDOS), K⁻ only on w_static.
"""
function fill_kp_lin!(Kp::Matrix{ComplexF64}, w, wp::Vector{Float64}, lk::LinearKernel)
    size(Kp) == (length(w), length(wp)) || throw(DimensionMismatch("Kp block does not match the grids"))
    @inbounds for i in eachindex(wp)
        b = wp[i]
        for j in eachindex(w)
            kp, _ = lin_kernel_eval(lk, float(w[j]), b)
            Kp[j, i] = kp
        end
    end
    return Kp
end

function fill_km_lin!(Km::Matrix{ComplexF64}, w, wp::Vector{Float64}, lk::LinearKernel)
    size(Km) == (length(w), length(wp)) || throw(DimensionMismatch("Km block does not match the grids"))
    @inbounds for i in eachindex(wp)
        b = wp[i]
        for j in eachindex(w)
            _, km = lin_kernel_eval(lk, float(w[j]), b)
            Km[j, i] = km
        end
    end
    return Km
end


"""
    lin_tail_accumulate_chi!(Ichi, A_dif, B_dif, A_sum, B_sum, noff, nout, wgt, fp, fm, g_chi)

K⁺-only tail accumulation for the χ channel (Ichi += ∫K⁺ g_chi over the tail). `noff` is the
difference-index offset of the shared tail arrays (their build length = master-grid length),
`nout` the number of output ω-points (here nout = noff = nchi).
"""
function lin_tail_accumulate_chi!(Ichi::Vector{ComplexF64},
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
    lin_tail_accumulate!(Iz, Iphi, A_dif, B_dif, A_sum, B_sum, ns, wgt, fp, fm, g_z, g_phi)

Accumulate the tail ω'-segment into Iz/Iphi from the precomputed samples — one loop over ω'
building the K_tail column for all ω by pure lookups (𝒦 = A + f·B), no interpolation.
`A_dif`/`B_dif` are indexed by p = (j-i)+noff, `A_sum`/`B_sum` by q = i+j-1. `noff` is the
difference offset of the shared tail arrays (their build/master length); `nout` is the number
of ω-output points (ns for the Km/Kp channels, which loop only the first nout ≤ noff rows).
"""
function lin_tail_accumulate!(Iz::Vector{ComplexF64}, Iphi::Vector{ComplexF64},
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
