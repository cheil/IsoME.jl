"""
    File containing the head/tail split of the ω'-integration grid for the
    real-axis vDOS solver.

    The grid is split at wp_max into
        - head: pole-anchored Chebyshev clusters in (0, wp_max], rebuilt whenever
                the integrand poles move (requires re-evaluating the kernel columns)
        - tail: fixed linear grid [wp_max, reOmega_c_shift], kernels evaluated
                once per temperature

    The ω'-integration itself is done as a matrix-vector product with the
    trapezoidal weights folded into the vector, so the (ω × ω') integrand
    matrices are never materialized.
"""


# pole movement (as fraction of wp_max) that triggers a rebuild of the head grid
const WPRIME_REFRESH_TOL = 1e-3


"""
    trapz_weights(x)

Trapezoidal quadrature weights on the grid x, such that dot(w, y) == trapz(x, y).
"""
function trapz_weights(x::AbstractVector{<:Real})
    n = length(x)
    wgt = zeros(n)
    @inbounds for i in 1:n-1
        h = (x[i+1] - x[i]) / 2
        wgt[i] += h
        wgt[i+1] += h
    end
    return wgt
end


"""
    cheb_halfgrid(m)

m Chebyshev-type nodes in (0, 1), ascending, clustered towards 0.
Same construction as the upper half of the grids in make_vDOS_wprime_grid.
"""
function cheb_halfgrid(m::Int)
    j = m:(2*m-1)
    x = 1 .+ cos.((2 .* j .+ 1) .* π ./ (4 * m))
    return reverse(x)
end


"""
    make_head_grid(poles, wp_max, pts_per_pole)

Head part of the ω'-grid: one Chebyshev cluster of `pts_per_pole` points per pole,
dense towards the pole from both sides and spanning up to the midpoints between
neighbouring poles (0 and wp_max at the outer ends). The clusters tile the whole
(0, wp_max] interval, so no separate background grid is needed.

Returns `K * pts_per_pole` strictly increasing points (K = number of distinct
poles inside the interval), with `pts[end] == wp_max` (the junction with the tail
grid). The length therefore varies with the number of poles.
"""
function make_head_grid(poles::Vector{Float64}, wp_max::Float64, pts_per_pole::Int)
    pts_per_pole >= 2 || error("pts_per_pole must be >= 2 (got $pts_per_pole)")

    # keep only poles strictly inside (0, wp_max), merge near-duplicates
    p_sorted = sort(filter(x -> isfinite(x) && 0 < x < wp_max, poles))
    tol = max(1e-3 * wp_max, 5) # poles min 5 meV apart
    p = [p_sorted[i] for i in eachindex(p_sorted) if i == 1 || p_sorted[i] - p_sorted[i-1] > tol]
    if isempty(p)
        p = [wp_max / 4]
    end
    K = length(p)

    # cluster boundaries: midpoints between neighbouring poles, 0 and wp_max outside
    bounds = [0.0; (p[1:end-1] .+ p[2:end]) ./ 2; wp_max]

    pts = Vector{Float64}()
    sizehint!(pts, K * pts_per_pole)
    for k in 1:K
        m_lo = pts_per_pole ÷ 2
        m_hi = pts_per_pole - m_lo
        h_lo = p[k] - bounds[k]
        h_hi = bounds[k+1] - p[k]
        append!(pts, p[k] .- h_lo .* reverse(cheb_halfgrid(m_lo)))   # (bounds[k], p), dense at p
        append!(pts, p[k] .+ h_hi .* cheb_halfgrid(m_hi))            # (p, bounds[k+1]), dense at p
    end

    sort!(pts)
    # enforce strictly increasing points (adjacent clusters touch at the boundaries)
    eps_min = 1e-12 * wp_max
    @inbounds for i in 2:length(pts)
        if pts[i] - pts[i-1] < eps_min
            pts[i] = pts[i-1] + eps_min
        end
    end
    pts[end] = wp_max   # junction with the tail grid must be exact

    @assert length(pts) == K * pts_per_pole
    return pts
end


"""
    grid_sign_roots(x, y)

Roots of y(x) located by sign changes on the grid, refined by linear interpolation.
"""
function grid_sign_roots(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    roots = Float64[]
    @inbounds for i in 1:length(y)-1
        y1, y2 = y[i], y[i+1]
        if isfinite(y1) && isfinite(y2) && ((y1 < 0) != (y2 < 0)) && y1 != y2
            push!(roots, x[i] - y1 * (x[i+1] - x[i]) / (y2 - y1))
        end
    end
    return roots
end


"""
    find_integrand_poles(wp, w_static, w_static_chi, znormip, deltaip, shiftip, fermi_level, gap0)

Locate the poles of the ε-integrated ω'-integrand from the previous iteration's
self-energy components:
    - roots of Im(χ ± ε_p)  (collapsing Lorentzian widths)
    - minimum of |ε_p|²     (branch point / gap edge)
Returns a sorted, non-empty vector (falls back to gap0 if nothing is found).
"""
function find_integrand_poles(wp::Vector{Float64}, w_static, w_static_chi,
                              znormip::Vector{ComplexF64}, deltaip::Vector{ComplexF64},
                              shiftip::Vector{ComplexF64}, fermi_level::Float64, gap0::Float64)

    Z_ongrid = linear_interpolation(w_static, znormip, extrapolation_bc=Flat()).(wp)
    phi_ongrid = linear_interpolation(w_static, deltaip .* znormip, extrapolation_bc=Flat()).(wp)
    shift_ongrid = linear_interpolation(w_static_chi, shiftip, extrapolation_bc=Flat()).(wp) .- fermi_level

    ε_p = sqrt.(wp .^ 2 .* Z_ongrid .^ 2 .- phi_ongrid .^ 2)

    poles = Float64[]
    #append!(poles, grid_sign_roots(wp, imag.(shift_ongrid .+ ε_p)))
    #append!(poles, grid_sign_roots(wp, imag.(shift_ongrid .- ε_p)))
    push!(poles, wp[argmin(abs2.(ε_p))])
    #print(poles)

    filter!(x -> isfinite(x) && x > 0, poles)
    isempty(poles) && push!(poles, gap0)
    sort!(poles)
    return poles
end


"""
    WPrimeWorkspace

ω'-grids, trapezoidal weights and precomputed kernel data for the head/tail split.
The junction point wp_max is contained in both head and tail, so
trapz(wp_full) == trapz(head) + trapz(tail) exactly.

The ω'-integral is done entirely through the 1-D linear kernel `lin` (𝒦 = A(x) + f(ω')·B(x)),
which needs only O(N) storage:
  * the materialized head blocks `Km_head_lin`/`Kp_head_lin` (ns×n_head resp. nmaster×n_head),
    so the pole-clustered head is a plain gemv;
  * the Fermi factors fp/fm = f(±ω') on the tail together with the Toeplitz/Hankel A,B samples
    (`A_dif_tail`, …), which let the uniform tail be summed on the fly (A + f·B) without ever
    materializing the O(N²) tail matrix.
`w_static` is a prefix of `w_static_chi` (equal step), so the K⁺ blocks live on the master
(χ) grid and serve both Iphi (first ns rows) and Ichi (all rows).
"""
mutable struct WPrimeWorkspace
    wp_max::Float64
    n_head::Int
    wp_head::Vector{Float64}
    wp_tail::Vector{Float64}
    wp_full::Vector{Float64}
    wgt_head::Vector{Float64}
    wgt_tail::Vector{Float64}
    Km_head_lin::Matrix{ComplexF64}     # K⁻(w_static, wp_head)   materialized head (gemv)
    Kp_head_lin::Matrix{ComplexF64}     # K⁺(w_master, wp_head)
    fp_tail::Vector{Float64}            # f(ω')  on wp_tail
    fm_tail::Vector{Float64}            # f(-ω') on wp_tail
    A_dif_tail::Vector{ComplexF64}      # A,B sampled on the tail difference/sum arguments
    B_dif_tail::Vector{ComplexF64}      # (Toeplitz/Hankel, O(ns + n_tail)); tail lookup route
    A_sum_tail::Vector{ComplexF64}
    B_sum_tail::Vector{ComplexF64}
    poles_prev::Vector{Float64}
    lin::LinearKernel                   # 1-D difference-grid kernel (A(x), B(x))
end

function WPrimeWorkspace(wp_head::Vector{Float64}, wp_tail::Vector{Float64},
                         w_static, w_static_chi;
                         poles::Vector{Float64}=[NaN], lin::LinearKernel)
    wp_head[end] == wp_tail[1] || error("head and tail grids must share the junction point")
    n_head = length(wp_head)
    ns = length(w_static)

    # master ω-grid for the K⁺ channel: w_static is a prefix of w_static_chi (equal step), so one
    # K⁺ block on the master grid covers Iphi (its first ns rows) and Ichi (all rows). K⁻ is only
    # needed on w_static. The tail A,B samples are likewise built once on the master grid.
    # w_static_chi === nothing -> cDOS mode: master grid == w_static.
    w_master = isnothing(w_static_chi) ? w_static : w_static_chi
    nmaster = length(w_master)

    fp_tail, fm_tail = fermi_pm(lin.β, wp_tail)
    A_dif_tail, B_dif_tail, A_sum_tail, B_sum_tail = precompute_tail_AB(lin, w_master, wp_tail)

    Km_head_lin = Matrix{ComplexF64}(undef, ns, n_head)       # K⁻ on w_static
    Kp_head_lin = Matrix{ComplexF64}(undef, nmaster, n_head)  # K⁺ on master grid
    fill_km_lin!(Km_head_lin, w_static, wp_head, lin)
    fill_kp_lin!(Kp_head_lin, w_master, wp_head, lin)

    return WPrimeWorkspace(wp_head[end], n_head, wp_head, wp_tail, vcat(wp_head, wp_tail),
                           trapz_weights(wp_head), trapz_weights(wp_tail),
                           Km_head_lin, Kp_head_lin,
                           fp_tail, fm_tail,
                           A_dif_tail, B_dif_tail, A_sum_tail, B_sum_tail,
                           copy(poles), lin)
end


"""
    lin_kernel_omega_integral(ws, w_static, g_z, g_phi)

O(N)-memory linear-route ω'-integrals ∫K⁻g_z and ∫K⁺g_phi over the full (head+tail) grid,
split into a head part and a tail part. The Chebyshev head uses the materialized
`Kp_head_lin`/`Km_head_lin` blocks (ns×n_head) via a gemv; the uniform tail uses the
precomputed Toeplitz/Hankel A,B samples (pure lookups) so it never materializes the O(N²)
tail matrix. Returns (Iz, Iphi).
"""
function lin_kernel_omega_integral(ws::WPrimeWorkspace, w_static, g_z, g_phi)
    nh = ws.n_head
    ntot = length(g_z)
    ns = length(w_static)
    noff = size(ws.Kp_head_lin, 1)      # master-grid length = difference offset of the tail arrays

    # head: materialized K blocks (gemv). K⁺ lives on the master grid; Iphi is its first ns rows.
    gh_z = ComplexF64.(ws.wgt_head .* view(g_z, 1:nh))
    gh_phi = ComplexF64.(ws.wgt_head .* view(g_phi, 1:nh))
    Iz = ws.Km_head_lin * gh_z
    Iphi = view(ws.Kp_head_lin, 1:ns, :) * gh_phi

    # tail: precomputed Toeplitz/Hankel A,B lookups (no matrix), loop only the first ns rows
    lin_tail_accumulate!(Iz, Iphi, ws.A_dif_tail, ws.B_dif_tail, ws.A_sum_tail, ws.B_sum_tail,
                         noff, ns, ws.wgt_tail, ws.fp_tail, ws.fm_tail,
                         view(g_z, nh+1:ntot), view(g_phi, nh+1:ntot))

    return Iz, Iphi
end


"""
    lin_kernel_omega_integral_chi(ws, w_static_chi, g_chi)

O(N) linear-route χ integral ∫K⁺ g_chi on the master (χ) grid: K⁺ head gemv on the full
`Kp_head_lin` plus the shared Toeplitz/Hankel tail. Returns Ichi (length nchi).
"""
function lin_kernel_omega_integral_chi(ws::WPrimeWorkspace, w_static_chi, g_chi)
    nh = ws.n_head
    ntot = length(g_chi)
    nchi = length(w_static_chi)
    noff = size(ws.Kp_head_lin, 1)      # = nchi

    gh_chi = ComplexF64.(ws.wgt_head .* view(g_chi, 1:nh))
    Ichi = ws.Kp_head_lin * gh_chi

    lin_tail_accumulate_chi!(Ichi, ws.A_dif_tail, ws.B_dif_tail, ws.A_sum_tail, ws.B_sum_tail,
                             noff, nchi, ws.wgt_tail, ws.fp_tail, ws.fm_tail,
                             view(g_chi, nh+1:ntot))

    return Ichi
end


"""
    build_wprime_workspace(inp, realAxisParameter, gap0_start)

Set up the vDOS workspace for one temperature. The head/tail split is at
wp_max = 2·Δ(0) of the starting state (the cDOS solution / previous temperature),
kept within [10·domega, reOmega_c_shift/2].
"""
function build_wprime_workspace(inp::arguments, realAxisParameter, gap0_start::Float64)
    (w_static, w_static_chi, lin_kernel) = realAxisParameter

    gap0 = (isfinite(gap0_start) && gap0_start > 0) ? gap0_start : inp.minGap
    wp_max = clamp(2 * gap0, 10 * inp.domega, inp.reOmega_c_shift / 2)
    pts_per_pole = inp.n_cheb    # Chebyshev points per pole in the head region

    # tail on the domega grid (same step as w_static and w_static_chi) -> all channels k=1
    wp_tail = collect(wp_max:inp.domega:inp.reOmega_c_shift)
    wp_head = make_head_grid([gap0], wp_max, pts_per_pole)

    return WPrimeWorkspace(wp_head, wp_tail, w_static, w_static_chi;
                           poles=[gap0], lin=lin_kernel)
end


"""
    find_integrand_poles_vDOSW(wp, w_static, w_static_chi, znormip, phi_mat, shiftip,
                               fermi_level, dos_en, idx_ef, gap0)

Poles of the ε-integrated ω'-integrand for the vDOS+W approximation, evaluated at the
Fermi-level ε-interval [dos_en[idx_ef], dos_en[idx_ef+1]]. Unlike the vDOS+μ case the
gap enters through the ε-dependent φ(ω,ε), so the pole positions follow the *modified*
quantities χ_mod = S and ε_p_mod = P of `epsilon_helpers_vDOSW`:

    S = (χ + Φ0·Φ1) / (1+Φ1²),   P = sqrt((ω²Z² − χ² − Φ0²)/(1+Φ1²) + S²)

The poles are the roots of Im(S ± P) (collapsing Lorentzian widths) plus the minimum of
|P|² (branch point). Returns a sorted, non-empty vector.
"""
function find_integrand_poles_vDOSW(wp::Vector{Float64}, w_static, w_static_chi,
                                    znormip::Vector{ComplexF64}, phi_mat::Matrix{ComplexF64},
                                    shiftip::Vector{ComplexF64}, fermi_level::Float64,
                                    dos_en::Vector{Float64}, idx_ef::Int, gap0::Float64)

    Z_ong = linear_interpolation(w_static, znormip, extrapolation_bc=Flat()).(wp)
    chi_ong = linear_interpolation(w_static_chi, shiftip, extrapolation_bc=Flat()).(wp) .- fermi_level

    # Fermi-level ε-interval and the affine coefficients of φ(ε) there
    j = clamp(idx_ef, 1, size(phi_mat, 2) - 1)
    dε = dos_en[j+1] - dos_en[j]
    phi_j = linear_interpolation(w_static, phi_mat[:, j], extrapolation_bc=Flat()).(wp)
    phi_jp = linear_interpolation(w_static, phi_mat[:, j+1], extrapolation_bc=Flat()).(wp)
    Φ1 = dε == 0 ? zero(phi_j) : (phi_jp .- phi_j) ./ dε   # keep ComplexF64 in both branches (else Φ1 is a type-unstable Union)
    Φ0 = phi_j .- dos_en[j] .* Φ1

    scale = 1 .+ Φ1 .^ 2
    wZ = wp .* Z_ong
    S = (chi_ong .+ Φ0 .* Φ1) ./ scale
    P = sqrt.((wZ .^ 2 .- chi_ong .^ 2 .- Φ0 .^ 2) ./ scale .+ S .^ 2)

    poles = Float64[]
    #append!(poles, grid_sign_roots(wp, imag.(S .+ P)))
    #append!(poles, grid_sign_roots(wp, imag.(S .- P)))
    push!(poles, wp[argmin(abs2.(P))])

    #println("Poles: ", poles)

    filter!(x -> isfinite(x) && x > 0, poles)
    isempty(poles) && push!(poles, gap0)
    sort!(poles)
    return poles
end


"""
    find_cDOS_pole(w_static, deltaip, gap0)

Pole of the cDOS+μ ω'-integrand: the root of Θ(ω') = ω' − Re Δ(ω') (the gap edge where
ω'² = Δ(ω')²). Falls back to gap0 if no root is found.
"""
function find_cDOS_pole(w_static, deltaip::Vector{ComplexF64}, gap0::Float64)
    Delta_func = linear_interpolation(w_static, deltaip, extrapolation_bc=Flat())
    root_eq = x -> x - real(Delta_func(x))
    for x0 in (gap0, 20.0)
        try
            return find_zero(root_eq, x0)
        catch
        end
    end
    #printWarning("Couldn't find the root of x-Δ(x). Shfiting chebyshev grid to BCS gap instead.")
    return gap0
end


"""
    build_cDOS_wprime_workspace(inp, realAxisParameter, gap0_start)

Head/tail ω'-workspace for the cDOS+μ approximation. The head/tail split is at
wp_max = 2·gap0_start (kept within [10·domega, reOmega_c/2]); the head is a Chebyshev
cluster around the pole of Θ(ω') and the tail is a static linear grid up to reOmega_c
(step domega) whose kernels are built once. Only the Km/Kp channels are built (no χ channel).
"""
function build_cDOS_wprime_workspace(inp::arguments, realAxisParameter, gap0_start::Float64)
    (w_static, _, lin_kernel) = realAxisParameter

    gap = (isfinite(gap0_start) && gap0_start > 0) ? gap0_start : inp.minGap
    wp_max = clamp(2 * gap, 10 * inp.domega, inp.reOmega_c / 2)
    pts_per_pole = inp.n_cheb    # Chebyshev points in the pole-anchored head region

    wp_tail = collect(wp_max:inp.domega:inp.reOmega_c)     # static linear tail, kernels built once
    wp_head = make_head_grid([gap], wp_max, pts_per_pole)  # Chebyshev cluster around the pole

    return WPrimeWorkspace(wp_head, wp_tail, w_static, nothing;
                           poles=[gap], lin=lin_kernel)
end


"""
    maybe_refresh_head!(ws, poles, w_static, w_static_chi, pts_per_pole)

Rebuild the head grid and its kernel columns if the poles moved by more than
WPRIME_REFRESH_TOL·wp_max since the last rebuild, or if their number changed.
Each pole gets `pts_per_pole` Chebyshev points, so the head length varies with the
number of poles; the head kernel buffers are reallocated when that length changes.
Returns true if rebuilt.
"""
function maybe_refresh_head!(ws::WPrimeWorkspace, poles::Vector{Float64},
                             w_static, w_static_chi, pts_per_pole::Int)
    stale = length(poles) != length(ws.poles_prev) ||
            any(!isfinite, ws.poles_prev) ||
            maximum(abs.(poles .- ws.poles_prev)) > WPRIME_REFRESH_TOL * ws.wp_max
    stale || return false

    #println("Refresh head: ", length(ws.poles_prev)," --> ",length(poles))

    ws.wp_head = make_head_grid(poles, ws.wp_max, pts_per_pole)
    ws.n_head = length(ws.wp_head)
    ws.wgt_head = trapz_weights(ws.wp_head)
    ws.wp_full = vcat(ws.wp_head, ws.wp_tail)

    # head length changes with the number of poles -> reallocate the head buffers
    # (row counts are preserved: K⁻ on w_static, K⁺ on the master grid)
    if size(ws.Km_head_lin, 2) != ws.n_head
        ws.Km_head_lin = Matrix{ComplexF64}(undef, size(ws.Km_head_lin, 1), ws.n_head)
        ws.Kp_head_lin = Matrix{ComplexF64}(undef, size(ws.Kp_head_lin, 1), ws.n_head)
    end

    @time w_master = isnothing(w_static_chi) ? w_static : w_static_chi
    fill_km_lin!(ws.Km_head_lin, w_static, ws.wp_head, ws.lin)   # K⁻ on w_static
    fill_kp_lin!(ws.Kp_head_lin, w_master, ws.wp_head, ws.lin)   # K⁺ on master grid
    ws.poles_prev = copy(poles)
    return true
end


"""
    wprime_trapz(ws, g)

Scalar trapezoidal integral of g over the full ω'-grid of the workspace.
"""
function wprime_trapz(ws::WPrimeWorkspace, g::AbstractVector)
    nh = ws.n_head
    return dot(ws.wgt_head, view(g, 1:nh)) + dot(ws.wgt_tail, view(g, nh+1:length(g)))
end


