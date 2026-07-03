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
    p = sort(filter(x -> isfinite(x) && 0 < x < wp_max, poles))
    tol = max(1e-3 * wp_max, 5) # poles min 5 meV apart
    p = [p[i] for i in eachindex(p) if i == 1 || p[i] - p[i-1] > tol]
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
    append!(poles, grid_sign_roots(wp, imag.(shift_ongrid .+ ε_p)))
    append!(poles, grid_sign_roots(wp, imag.(shift_ongrid .- ε_p)))
    push!(poles, wp[argmin(abs2.(ε_p))])

    filter!(x -> isfinite(x) && x > 0, poles)
    isempty(poles) && push!(poles, gap0)
    sort!(poles)
    return poles
end


"""
    evaluate_Kernels!(K_vals, A, B, K_func)

In-place, threaded version of evaluate_Kernels.
"""
function evaluate_Kernels!(K_vals::Matrix{ComplexF64}, A, B, K_func)
    size(K_vals) == (length(A), length(B)) || throw(DimensionMismatch("kernel buffer does not match the grids"))
    # Add parallelization 
    #Threads.@threads for i in eachindex(B)
    for i in eachindex(B)
        b = B[i]
        @inbounds for j in eachindex(A)
            K_vals[j, i] = K_func(A[j], b)
        end
    end
    return K_vals
end


"""
    WPrimeWorkspace

ω'-grids, trapezoidal weights and precomputed kernel blocks for the head/tail split.
The junction point wp_max is contained in both head and tail, so
trapz(wp_full) == trapz(head) + trapz(tail) exactly.
"""
mutable struct WPrimeWorkspace
    wp_max::Float64
    n_head::Int
    wp_head::Vector{Float64}
    wp_tail::Vector{Float64}
    wp_full::Vector{Float64}
    wgt_head::Vector{Float64}
    wgt_tail::Vector{Float64}
    Km_head::Matrix{ComplexF64}         # K⁻(w_static, wp_head)
    Kp_head::Matrix{ComplexF64}         # K⁺(w_static, wp_head)
    Kpchi_head::Matrix{ComplexF64}      # K⁺(w_static_chi, wp_head)
    Km_tail::Matrix{ComplexF64}
    Kp_tail::Matrix{ComplexF64}
    Kpchi_tail::Matrix{ComplexF64}
    poles_prev::Vector{Float64}
end

function WPrimeWorkspace(wp_head::Vector{Float64}, wp_tail::Vector{Float64},
                         w_static, w_static_chi, Kp_func, Km_func;
                         poles::Vector{Float64}=[NaN])
    wp_head[end] == wp_tail[1] || error("head and tail grids must share the junction point")
    n_head = length(wp_head)
    n_tail = length(wp_tail)
    # w_static_chi === nothing -> cDOS mode: no χ channel, empty (0-row) buffers
    nchi = isnothing(w_static_chi) ? 0 : length(w_static_chi)

    ws = WPrimeWorkspace(wp_head[end], n_head, wp_head, wp_tail, vcat(wp_head, wp_tail),
                         trapz_weights(wp_head), trapz_weights(wp_tail),
                         Matrix{ComplexF64}(undef, length(w_static), n_head),
                         Matrix{ComplexF64}(undef, length(w_static), n_head),
                         Matrix{ComplexF64}(undef, nchi, n_head),
                         Matrix{ComplexF64}(undef, length(w_static), n_tail),
                         Matrix{ComplexF64}(undef, length(w_static), n_tail),
                         Matrix{ComplexF64}(undef, nchi, n_tail),
                         copy(poles))

    evaluate_Kernels!(ws.Km_head, w_static, wp_head, Km_func)
    evaluate_Kernels!(ws.Kp_head, w_static, wp_head, Kp_func)
    evaluate_Kernels!(ws.Km_tail, w_static, wp_tail, Km_func)
    evaluate_Kernels!(ws.Kp_tail, w_static, wp_tail, Kp_func)
    if nchi > 0
        evaluate_Kernels!(ws.Kpchi_head, w_static_chi, wp_head, Kp_func)
        evaluate_Kernels!(ws.Kpchi_tail, w_static_chi, wp_tail, Kp_func)
    end

    return ws
end


"""
    build_wprime_workspace(inp, realAxisParameter, gap0_start)

Set up the vDOS workspace for one temperature. The head/tail split is at
wp_max = 2·Δ(0) of the starting state (the cDOS solution / previous temperature),
kept within [10·domega_shift, reOmega_c_shift/2].
"""
function build_wprime_workspace(inp::arguments, realAxisParameter, gap0_start::Float64)
    (Kp_func, Km_func, w_static, _, w_static_chi) = realAxisParameter

    gap0 = (isfinite(gap0_start) && gap0_start > 0) ? gap0_start : inp.minGap
    wp_max = clamp(2 * gap0, 10 * inp.domega_shift, inp.reOmega_c_shift / 2)
    pts_per_pole = inp.n_cheb    # Chebyshev points per pole in the head region

    wp_tail = collect(wp_max:inp.domega_shift:inp.reOmega_c_shift)
    wp_head = make_head_grid([gap0], wp_max, pts_per_pole)

    return WPrimeWorkspace(wp_head, wp_tail, w_static, w_static_chi, Kp_func, Km_func; poles=[gap0])
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
    Φ1 = dε == 0 ? zero(wp) : (phi_jp .- phi_j) ./ dε
    Φ0 = phi_j .- dos_en[j] .* Φ1

    scale = 1 .+ Φ1 .^ 2
    wZ = wp .* Z_ong
    S = (chi_ong .+ Φ0 .* Φ1) ./ scale
    P = sqrt.((wZ .^ 2 .- chi_ong .^ 2 .- Φ0 .^ 2) ./ scale .+ S .^ 2)

    poles = Float64[]
    append!(poles, grid_sign_roots(wp, imag.(S .+ P)))
    append!(poles, grid_sign_roots(wp, imag.(S .- P)))
    push!(poles, wp[argmin(abs2.(P))])

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
    build_cDOS_wprime_workspace(inp, realAxisParameter, BCS_gap)

Head/tail ω'-workspace for the cDOS+μ approximation. The head/tail split is at
wp_max = 2·BCS_gap, the tail is a linear grid up to reOmega_c (step domega), and only the
Km/Kp channels are built (no χ channel).
"""
function build_cDOS_wprime_workspace(inp::arguments, realAxisParameter, BCS_gap::Float64)
    (Kp_func, Km_func, w_static, _, _) = realAxisParameter

    gap = (isfinite(BCS_gap) && BCS_gap > 0) ? BCS_gap : inp.minGap
    gap = inp.reOmega_c
    wp_max = clamp(2 * gap, 10 * inp.domega, inp.reOmega_c / 2)
    wp_max = inp.reOmega_c
    pts_per_pole = inp.n_cheb

    print(gap)

    wp_tail = collect(wp_max:inp.domega:inp.reOmega_c)
    wp_head = make_head_grid([gap], wp_max, pts_per_pole)

    # # Chebyshev grid over whole range, center it at pole of Θ at each iteration
    #w_dynam = reverse(inp.reOmega_c .+ inp.reOmega_c .* cos.((2 .* (inp.n_cheb/2:inp.n_cheb-1) .+ 1) .* π ./ (2 * inp.n_cheb)))    # only positvie chebyshev nodes
    #wp_head = w_dynam
    # # currently only uses w_dynam
    # wp_head = w_dynam
    # wp_tail = [w_dynam[end]]

    return WPrimeWorkspace(wp_head, wp_tail, w_static, nothing, Kp_func, Km_func; poles=[gap])
end


"""
    maybe_refresh_head!(ws, poles, w_static, w_static_chi, Kp_func, Km_func, pts_per_pole)

Rebuild the head grid and its kernel columns if the poles moved by more than
WPRIME_REFRESH_TOL·wp_max since the last rebuild, or if their number changed.
Each pole gets `pts_per_pole` Chebyshev points, so the head length varies with the
number of poles; the head kernel buffers are reallocated when that length changes.
Returns true if rebuilt.
"""
function maybe_refresh_head!(ws::WPrimeWorkspace, poles::Vector{Float64},
                             w_static, w_static_chi, Kp_func, Km_func, pts_per_pole::Int)
    stale = length(poles) != length(ws.poles_prev) ||
            any(!isfinite, ws.poles_prev) ||
            maximum(abs.(poles .- ws.poles_prev)) > WPRIME_REFRESH_TOL * ws.wp_max
    stale || return false

    ws.wp_head = make_head_grid(poles, ws.wp_max, pts_per_pole)
    ws.n_head = length(ws.wp_head)
    ws.wgt_head = trapz_weights(ws.wp_head)
    ws.wp_full = vcat(ws.wp_head, ws.wp_tail)

    # head length changes with the number of poles -> reallocate the head buffers
    # (row counts are preserved: 0 rows for the χ block keeps the cDOS mode χ-less)
    if size(ws.Km_head, 2) != ws.n_head
        ws.Km_head = Matrix{ComplexF64}(undef, size(ws.Km_head, 1), ws.n_head)
        ws.Kp_head = Matrix{ComplexF64}(undef, size(ws.Kp_head, 1), ws.n_head)
        ws.Kpchi_head = Matrix{ComplexF64}(undef, size(ws.Kpchi_head, 1), ws.n_head)
    end

    evaluate_Kernels!(ws.Km_head, w_static, ws.wp_head, Km_func)
    evaluate_Kernels!(ws.Kp_head, w_static, ws.wp_head, Kp_func)
    if size(ws.Kpchi_head, 1) > 0
        evaluate_Kernels!(ws.Kpchi_head, w_static_chi, ws.wp_head, Kp_func)
    end
    ws.poles_prev = copy(poles)
    return true
end


"""
    kernel_omega_integral(K_head, K_tail, ws, g)

∫ K(ω,ω') g(ω') dω' over the full ω'-grid as two gemv calls with the trapezoidal
weights folded into g. g is indexed like ws.wp_full.
"""
function kernel_omega_integral(K_head::Matrix{ComplexF64}, K_tail::Matrix{ComplexF64},
                               ws::WPrimeWorkspace, g::AbstractVector{<:Real})
    nh = ws.n_head
    gh = ComplexF64.(ws.wgt_head .* view(g, 1:nh))
    gt = ComplexF64.(ws.wgt_tail .* view(g, nh+1:length(g)))
    return K_head * gh .+ K_tail * gt
end


"""
    wprime_trapz(ws, g)

Scalar trapezoidal integral of g over the full ω'-grid of the workspace.
"""
function wprime_trapz(ws::WPrimeWorkspace, g::AbstractVector)
    nh = ws.n_head
    return dot(ws.wgt_head, view(g, 1:nh)) + dot(ws.wgt_tail, view(g, nh+1:length(g)))
end


##################################################################
# ---------------------- debug plotting ------------------------ #
##################################################################

# counter shown in the plot titles so successive iterations can be told apart
# (the files themselves are overwritten on every call)
const WPRIME_DEBUG_COUNT = Ref(0)

"""
    debug_plot_wprime(ws, integrands; outdir="")

Debug plots for the pole-anchored ω'-grid, enabled in realEliashbergEq via
ENV["ISOME_DEBUG_POLES"] = "1". Written to outdir (default: current directory):

    - wprime_debug_grid.png:  local grid spacing Δω' over the full ω'-grid
                              (head+tail), poles and wp_max marked
    - wprime_debug_<z|phi|chi>.png: integrand over the head region (top panel)
                              and local grid spacing (bottom panel), same x-axis,
                              to check that the grid clusters where the
                              integrands are peaked
    - wprime_debug_pole_<k>.png: zoom of all three integrands (normalized to
                              max |g| = 1 in the window) onto the grid cluster
                              around pole k

integrands is indexed like ws.wp_full (as inside realEliashbergEq).
"""
function debug_plot_wprime(ws::WPrimeWorkspace, integrands; outdir::String="")
    WPRIME_DEBUG_COUNT[] += 1
    it = WPRIME_DEBUG_COUNT[]
    poles = filter(isfinite, ws.poles_prev)
    nh = ws.n_head
    labels = ("z", "phi", "chi")

    # ---- overview: local spacing of the full grid ----
    pov = plot(ws.wp_head[2:end], diff(ws.wp_head); seriestype=:scatter, ms=1.5, msw=0,
               xscale=:log10, yscale=:log10, label="head",
               xlabel="ω' / meV", ylabel="local spacing Δω' / meV",
               title="ω'-grid spacing (call $it)", legend=:bottomright)
    plot!(pov, ws.wp_tail[2:end], diff(ws.wp_tail); seriestype=:scatter, ms=1.5, msw=0, label="tail")
    isempty(poles) || vline!(pov, poles; ls=:dash, lc=:red, label="poles")
    vline!(pov, [ws.wp_max]; ls=:dot, lc=:black, label="wp_max")
    savefig(pov, joinpath(outdir, "wprime_debug_grid.png"))

    # ---- per-integrand: integrand vs grid density in the head region ----
    for (idx, g) in enumerate(integrands)
        gh = real.(g[1:nh])

        ptop = plot(ws.wp_head, gh; label="integrand", ylabel="integrand ($(labels[idx]))",
                    title="ω'-integrand vs grid density (call $it)", xlims=(0, ws.wp_max))
        scatter!(ptop, ws.wp_head, gh; ms=1.2, msw=0, label="grid points")
        isempty(poles) || vline!(ptop, poles; ls=:dash, lc=:red, label="poles")

        pbot = plot(ws.wp_head[2:end], diff(ws.wp_head); seriestype=:scatter, ms=1.5, msw=0,
                    yscale=:log10, label="", xlabel="ω' / meV", ylabel="Δω' / meV",
                    xlims=(0, ws.wp_max))
        isempty(poles) || vline!(pbot, poles; ls=:dash, lc=:red, label="")

        savefig(plot(ptop, pbot; layout=(2, 1), link=:x),
                joinpath(outdir, "wprime_debug_$(labels[idx]).png"))
    end

    # ---- zoom onto the cluster around each pole (all integrands, normalized) ----
    for (k, p) in enumerate(poles)
        i0 = argmin(abs.(ws.wp_head .- p))
        window = max(1, i0 - 150):min(nh, i0 + 150)
        x = ws.wp_head[window]

        pz = plot(xlabel="ω' / meV", ylabel="integrand / max|integrand|",
                  title="pole $k at ω' = $(round(p, digits=4)) meV (call $it)")
        for (idx, g) in enumerate(integrands)
            gw = real.(g[window])
            gmax = maximum(abs, gw)
            gmax > 0 && plot!(pz, x, gw ./ gmax; label=labels[idx])
        end
        scatter!(pz, x, zero(x); ms=1.2, msw=0, label="grid points")
        vline!(pz, [p]; ls=:dash, lc=:red, label="pole")
        savefig(pz, joinpath(outdir, "wprime_debug_pole_$k.png"))
    end

    return nothing
end
