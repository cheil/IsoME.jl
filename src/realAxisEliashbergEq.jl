"""
    file containing the real axis eliashberg equations
"""

# --- DEBUG: per-iteration memory report (remove later) -----------------------
# report_iteration_mem(i_it, gridws; name1 = var1, ...) prints, once per iteration:
#   * the process peak RSS so far,
#   * the tracked variables ranked by live size, biggest first - the one that
#     dominates the working set is named in the header,
#   * a field-by-field breakdown of the ω'-workspace, which is usually that biggest
#     one. Only array-valued fields are listed (scalars cannot contribute
#     meaningfully); arrays inside nested struct fields are recursed into and shown
#     dotted, e.g. ker.A.
#
# The per-variable figures are separate Base.summarysize calls, so structure SHARED
# between two tracked variables is counted in both - realAxisParameter and gridws
# both reference the same PhononKernel, and that kernel shows up in full under each.
# `tracked` in the header is one summarysize over all of them at once and therefore
# does not double count; it is the figure to compare against peak RSS.
function report_iteration_mem(i_it, gridws; vars...)
    # top-level tracked variables, ranked by size (gridws is always tracked)
    tracked = Tuple{String,Int}[("gridws", Base.summarysize(gridws))]
    objs = Any[gridws]
    for (name, v) in vars
        push!(tracked, (string(name), Base.summarysize(v)))
        push!(objs, v)
    end
    sort!(tracked, by = t -> -t[2])
    unique_bytes = Base.summarysize(objs)

    # ω'-workspace field breakdown
    entries = Tuple{String,Int,String}[]     # (name, bytes, "dims eltype")
    collect_array_fields!(entries, gridws, "")
    sort!(entries, by = e -> -e[2])
    total = isempty(entries) ? 0 : sum(e -> e[2], entries)

    @printf("[MEM] it %-4d  peak RSS %.1f MiB   tracked %.1f MiB   biggest %s %.2f MiB\n",
            i_it, Sys.maxrss() / 2^20, unique_bytes / 2^20,
            tracked[1][1], tracked[1][2] / 2^20)
    for (name, b) in tracked
        @printf("      var %-22s %8.2f MiB\n", name, b / 2^20)
    end
    @printf("      gridws arrays %.2f MiB\n", total / 2^20)
    for (name, b, dims) in entries
        @printf("        %-14s %-26s %8.2f MiB  (%4.1f %%)\n",
                name, dims, b / 2^20, total == 0 ? 0.0 : 100 * b / total)
    end
    return nothing
end

# Walk `obj`s fields, collecting every AbstractArray; recurse into struct-valued
# fields so nested arrays (ker.A, ker.B) are reported separately rather than as
# one opaque lump. Numbers are skipped.
function collect_array_fields!(entries, obj, prefix)
    for name in fieldnames(typeof(obj))
        v = getfield(obj, name)
        v isa Number && continue
        label = prefix * string(name)
        if v isa AbstractArray
            push!(entries, (label, Base.summarysize(v),
                            join(size(v), "×") * " " * string(eltype(v))))
        elseif Base.isstructtype(typeof(v))
            collect_array_fields!(entries, v, label * ".")
        end
    end
    return entries
end


"""
    realEliashbergEq(beta, znormip, phi_ph_ip, phi_c_ip, shiftip, ws, w_static, dosef,
                     epsilon, dos, Weep, idx_ef, fermi_level, wgCoulomb, electronic_spec)

Real axis Eliashberg equations in the vDOS+W approximation
"""
function realEliashbergEq(beta::Float64, znormip::Vector{ComplexF64}, phi_ph_ip::Vector{ComplexF64}, phi_c_ip::Vector{ComplexF64}, shiftip::Vector{ComplexF64},
                          ws::WPrimeWorkspace, w_static::AbstractVector,
                          dosef::Float64, epsilon::Vector{Float64}, dos::Vector{Float64}, Weep::Matrix{Float64},
                          idx_ef::Int64, fermi_level::Float64, wgCoulomb::Float64, electronic_spec::Tuple)

    w_prime = ws.wp_full

    # interpolate Z,χ,ϕ onto ω'-integration grid
    Z_itp = linear_interpolation(w_static, znormip, extrapolation_bc=Flat())
    shift_itp = linear_interpolation(w_static, shiftip, extrapolation_bc=Flat())

    Z_ongrid = Z_itp.(w_prime)
    shift_ongrid = shift_itp.(w_prime) .- fermi_level
    # φ(ω,ε) = φ_ph(ω) + φ_c(ε): only φ_ph needs the ω→ω' interpolation; φ_c lives on the ε-grid.
    phi_ph_ongrid = linear_interpolation(w_static, phi_ph_ip, extrapolation_bc=Flat()).(w_prime)

    # ------------- ε-integration ------------- #
    integrands, coulomb_spectral = epsilon_helpers_vDOSW(electronic_spec, Weep, Z_ongrid, phi_ph_ongrid, phi_c_ip, shift_ongrid, w_prime)

    # DEBUG: the ε-integration allocates two M×Nint transients (C0phi/C1phi, freed on return)
    # plus the returned ndos×M coulomb_spectral; Sys.maxrss() is monotonic so it captures that
    # peak even though C0phi/C1phi are already gone by iteration-end. (remove later)
    # @info "[MEM] ε-integration" peak_rss_MiB = round(Sys.maxrss() / 2^20, digits = 1) coulomb_spectral_MiB = round(Base.summarysize(coulomb_spectral) / 2^20, digits = 1) C0C1phi_transient_MiB = round(2 * length(w_prime) * length(electronic_spec[5]) * 8 / 2^20, digits = 1)

    # Z-integrand must be positive (causality), flip the others accordingly
    # Note: Minus sign because split of -Θ to (ε+χ +- ε_p) (eq. (51))
    FLIP = ifelse.(integrands[1] .<= 0, 1, -1)
    g_z = -abs.(integrands[1])
    g_phi = -(FLIP .* integrands[2])
    g_chi = -(FLIP .* integrands[3])

    # ------------- Ω & ω'-integration ------------- #
    # O(N)-memory kernel route, one pass: K⁻ g_z plus K⁺ g_phi and K⁺ g_chi, all on w_static.
    Iz, Iphi, Ichi = kernel_omega_integral(ws, w_static, g_z, g_phi, g_chi)

    # Coulomb term: integrate N(ε')W(ε,ε') with the same piecewise-linear spectral
    # quadrature (coulomb_spectral is ndos × M). The dosef prefactor cancels.
    wgt_full = vcat(ws.wgt_head, ws.wgt_tail)
    coulomb_weight = wgt_full .* FLIP .* tanh.(beta .* w_prime ./ 2)
    phi_c_val = (wgCoulomb / π) .* (coulomb_spectral * coulomb_weight)   # φ_c(ε), length ndos

    Zval = 1 .+ Iz ./ (w_static .* (π * dosef))
    phi_ph_val = Iphi ./ (π * dosef)                                     # φ_ph(ω), length nω
    shift_val = -Ichi ./ (π * dosef)

    # φ(ω,ε) = φ_ph(ω) + φ_c(ε) is kept in separable form; no nω×ndos matrix is built.
    # Δ is derived from φ and Z in the solver after mixing; hand back the two φ vectors here.
    return Zval, shift_val, phi_ph_val, phi_c_val
end


"""
    realEliashbergEq(mu_star, beta, znormip, phiip, shiftip, ws, w_static, dosef, epsilon, dos, fermi_level)

Real axis Eliashberg equations in the vDOS+μ approximation. 
"""
function realEliashbergEq(mu_star::Float64, beta::Float64, znormip::Vector{ComplexF64}, phiip::Vector{ComplexF64},
                          shiftip::Vector{ComplexF64}, ws::WPrimeWorkspace, w_static::AbstractVector,
                          dosef::Float64, wgCoulomb::Float64,
                          fermi_level::Float64, electronic_spec::Tuple)

    w_prime = ws.wp_full

    # interpolate Z,χ,ϕ onto ω'-integration grid (φ is handed in directly)
    Z_itp = linear_interpolation(w_static, znormip, extrapolation_bc=Flat())
    phi_itp = linear_interpolation(w_static, phiip, extrapolation_bc=Flat())
    shift_itp = linear_interpolation(w_static, shiftip, extrapolation_bc=Flat())

    Z_ongrid = Z_itp.(w_prime)
    phi_ongrid = phi_itp.(w_prime)
    shift_ongrid = shift_itp.(w_prime) .- fermi_level

    # ------------- ε-integration ------------- #
    integrands = epsilon_helpers(electronic_spec, Z_ongrid, phi_ongrid, shift_ongrid, w_prime)

    # Z-integrand must be positive (causality), flip the others accordingly
    # Note: Minus sign because split of -Θ to (ε+χ +- ε_p) (eq. (51))
    FLIP = ifelse.(integrands[1] .<= 0, 1, -1)
    g_z = -abs.(integrands[1])
    g_phi = -(FLIP .* integrands[2])
    g_chi = -(FLIP .* integrands[3])

    # ------------- ω'-integration ------------- #
    # O(N)-memory kernel route, one pass: K⁻ g_z plus K⁺ g_phi and K⁺ g_chi, all on w_static.
    Iz, Iphi, Ichi = kernel_omega_integral(ws, w_static, g_z, g_phi, g_chi)

    # μ* Coulomb term
    coulomb = wprime_trapz(ws, g_phi .* tanh.(beta .* w_prime ./ 2))

    Zval = 1 .+ Iz ./ (w_static .* (π * dosef))
    phi_val = (Iphi .- wgCoulomb .* mu_star .* coulomb) ./ (π * dosef)
    shift_val = -Ichi ./ (π * dosef)

    return Zval, shift_val, phi_val
end




"""
     realEliashbergEq(mu_star, beta, deltaip, ws, w_static)

Real axis Eliashberg equations in cDOS+μ approximation
"""
function realEliashbergEq(mu_star::Float64, beta::Float64, deltaip::Vector{ComplexF64},
                          ws::WPrimeWorkspace, w_static::AbstractVector, wgCoulomb::Float64)

    Delta_func = linear_interpolation(w_static, deltaip, extrapolation_bc=Flat())

    w_prime = ws.wp_full
    Delta_ongrid = Delta_func.(w_prime)
    sqrt_eval = @. sqrt(w_prime^2 - Delta_ongrid^2)

    g_z = @. real(w_prime / sqrt_eval)          # Z   integrand density (∝ quasiparticle DOS)
    g_phi = @. real(Delta_ongrid / sqrt_eval)   # Δ·Z integrand density

    # O(N)-memory kernel route: materialized K head (gemv) + Toeplitz/Hankel tail (A + f·B on the fly)
    Iz, Iphi = kernel_omega_integral(ws, w_static, g_z, g_phi)

    # μ* Coulomb term (scalar, same for all ω)
    coulomb = wprime_trapz(ws, g_phi .* tanh.(beta .* w_prime ./ 2))

    Zval = 1 .- Iz ./ w_static
    Delta = (Iphi .- wgCoulomb .* mu_star .* coulomb) ./ Zval

    return Zval, Delta
end



###################################################
# ------------------- Helpers ------------------- #
###################################################
"""
    epsilon_helpers(electronic_spec, Z_ongrid, phi_ongrid, shift_ongrid, w_prime; idx_skip = nothing)

Sovle the epsilon integral in the vDOS+μ approximation
Returns the ω'-integrands
"""
function epsilon_helpers(electronic_spec, Z_ongrid, phi_ongrid, shift_ongrid, w_prime; idx_skip = nothing)
    eps_1, eps_2, dos_1, dos_2, dε, ddos, _, _ = electronic_spec

    # assuming a linear form of the dos
    M1 = transpose((dos_2 - dos_1)./(eps_2 - eps_1))
    M0 = transpose(dos_1) .- transpose(eps_1) .* M1
    
    # Helper functions for the ε-integration
    ε_p = sqrt.(w_prime.^2 .*Z_ongrid.^2 .- phi_ongrid.^2)       
    Rplus = real.(shift_ongrid .+ ε_p)
    Rminus = real.(shift_ongrid .- ε_p)
    Iplus = imag.(shift_ongrid .+ ε_p)
    Iminus = imag.(shift_ongrid .- ε_p)

    # M0- and M1-weighted ε-sums of the Lorentzian moments (pp and mm branches)
    lorentzian_moments_p = lorentzian_moments_weighted(eps_1, eps_2, Rplus, Iplus, M0, M1)
    lorentzian_moments_m = lorentzian_moments_weighted(eps_1, eps_2, Rminus, Iminus, M0, M1)

    # Evaluate ε-integrals to get ω'-integrands 
    integrands = Vector{Vector{Float64}}(undef, 3)
    for (idx, g) in enumerate((w_prime .* Z_ongrid, phi_ongrid, shift_ongrid))
        if ~isnothing(idx_skip) && any(idx .== idx_skip)
            continue
        end
        
        integrands[idx] = eval_spectral_integrals(g, lorentzian_moments_p, lorentzian_moments_m, ε_p, Rplus, Rminus, Iplus, Iminus, idx == 3)
    end

    return integrands
end

"""
    lorentzian_interval_moments(εl, εr, a, b)

Scalar Lorentzian moments ∫ εᵏ/((ε+a)²+b²) dε over [εl, εr] for k = 0..3.
for a single (ω′, ε-interval) pair.
"""
@inline function lorentzian_interval_moments(εl, εr, a, b)
    b2 = b^2
    a2 = a^2

    i0l = atan((a + εl) / b) / b
    i0r = atan((a + εr) / b) / b

    h_l = log((a + εl)^2 + b2)
    h_r = log((a + εr)^2 + b2)

    i1l = 0.5 * h_l - a * i0l
    i1r = 0.5 * h_r - a * i0r

    i2l = -a * h_l + (a2 - b2) * i0l + εl
    i2r = -a * h_r + (a2 - b2) * i0r + εr

    i3l = i0l * (3 * a * b2 - a * a2) + 0.5 * (3 * a2 - b2) * h_l + 0.5 * εl * (εl - 4 * a)
    i3r = i0r * (3 * a * b2 - a * a2) + 0.5 * (3 * a2 - b2) * h_r + 0.5 * εr * (εr - 4 * a)

    return i0r - i0l, i1r - i1l, i2r - i2l, i3r - i3l
end

"""
    spectral_C0C1(Rg, Ig, Rp, Rm, Ip, Im_, i0pp, i1pp, i2pp, i0mm, i1mm, i2mm)

Scalar version of `eval_spectral_integral_coefficients` for one (ω′, ε-interval) pair.
The slope (g_slope) contribution has the same functional form with the moments shifted
by one order, so it is obtained by calling this with (i1, i2, i3) instead of (i0, i1, i2).
"""
@inline function spectral_C0C1(Rg, Ig, Rp, Rm, Ip, Im_, i0pp, i1pp, i2pp, i0mm, i1mm, i2mm)
    C0 = Ig * (i1pp - i1mm) + Ig * (Rp * i0pp - Rm * i0mm) + Rg * (-Ip * i0pp + Im_ * i0mm)
    C1 = Ig * (i2pp - i2mm) + Ig * (Rp * i1pp - Rm * i1mm) + Rg * (-Ip * i1pp + Im_ * i1mm)
    return C0, C1
end

"""
    epsilon_helpers_vDOSW(electronic_spec, Weep, Z_ongrid, phi_ph_ongrid, phi_c, shift_ongrid, w_prime)

Fused route for the vDOS+W spectral integrals: instead of materialising the
I⁰..I³ Lorentzian-moment matrices (and S, P, R±, I± as full ω′×ε matrices), 
loop over the ε-intervals and accumulate the M0/M1-weighted ε-sum directly 
into length-M vectors for the Z, φ and χ integrands.
Only the φ-branch coefficients C0/C1 are stored as matrices, because the Coulomb
part couples them to M0W/M1W over the *output* ε-grid via a GEMM.
"""
function epsilon_helpers_vDOSW(electronic_spec, Weep, Z_ongrid, phi_ph_ongrid, phi_c, shift_ongrid, w_prime)
    eps_1, eps_2, dos_1, dos_2, dε, ddos, M0W, M1W = electronic_spec

    M = length(w_prime)
    Nint = length(dε)

    z_int = zeros(Float64, M)
    phi_int = zeros(Float64, M)
    chi_int = zeros(Float64, M)

    # φ-branch coefficients: needed as full matrices for the Coulomb GEMM below
    C0phi = Matrix{Float64}(undef, M, Nint)
    C1phi = Matrix{Float64}(undef, M, Nint)

    wZ = w_prime .* Z_ongrid

    @inbounds for jε in 1:Nint
        xl = eps_1[jε]
        xr = eps_2[jε]
        invdε = 1.0 / dε[jε]
        m1 = ddos[jε] * invdε
        m0 = dos_1[jε] - xl * m1

        # φ(ω,ε) = φ_ph(ω) + φ_c(ε): the ε-slope Φ1 = φ_c'(ε) is ω-independent
        Φ1 = (phi_c[jε+1] - phi_c[jε]) * invdε
        inv_scale = 1 / (1 + Φ1^2)

        for iw in 1:M
            Φ0 = phi_ph_ongrid[iw] + phi_c[jε] - xl * Φ1
            sh = shift_ongrid[iw]
            wz = wZ[iw]

            S = (sh + Φ0 * Φ1) * inv_scale
            P = sqrt((wz^2 - sh^2 - Φ0^2) * inv_scale + S^2)

            Sp = S + P
            Sm = S - P
            Rp = real(Sp); Ip = imag(Sp)
            Rm = real(Sm); Im_ = imag(Sm)

            i0pp, i1pp, i2pp, i3pp = lorentzian_interval_moments(xl, xr, Rp, Ip)
            i0mm, i1mm, i2mm, i3mm = lorentzian_interval_moments(xl, xr, Rm, Im_)

            c = inv_scale / (2 * P)

            # Z branch: g = ωZ/scale, no slope
            gz = wz * c
            C0, C1 = spectral_C0C1(real(gz), imag(gz), Rp, Rm, Ip, Im_,
                                   i0pp, i1pp, i2pp, i0mm, i1mm, i2mm)
            z_int[iw] += m0 * C0 + m1 * C1

            # φ branch: g = Φ0/scale, slope = Φ1/scale
            gphi = Φ0 * c
            gs = Φ1 * c
            C0, C1 = spectral_C0C1(real(gphi), imag(gphi), Rp, Rm, Ip, Im_,
                                   i0pp, i1pp, i2pp, i0mm, i1mm, i2mm)
            D0, D1 = spectral_C0C1(real(gs), imag(gs), Rp, Rm, Ip, Im_,
                                   i1pp, i2pp, i3pp, i1mm, i2mm, i3mm)
            C0 += D0
            C1 += D1
            C0phi[iw, jε] = C0
            C1phi[iw, jε] = C1
            phi_int[iw] += m0 * C0 + m1 * C1

            # χ branch: g = χ/scale, slope = 1/scale
            gchi = sh * c
            C0, C1 = spectral_C0C1(real(gchi), imag(gchi), Rp, Rm, Ip, Im_,
                                   i0pp, i1pp, i2pp, i0mm, i1mm, i2mm)
            D0, D1 = spectral_C0C1(real(c), imag(c), Rp, Rm, Ip, Im_,
                                   i1pp, i2pp, i3pp, i1mm, i2mm, i3mm)
            chi_int[iw] += m0 * (C0 + D0) + m1 * (C1 + D1)
        end
    end

    coulomb_spectral = Matrix{Float64}(undef, size(M0W, 1), M)
    mul!(coulomb_spectral, M0W, transpose(C0phi))
    mul!(coulomb_spectral, M1W, transpose(C1phi), 1.0, 1.0)

    return [z_int, phi_int, chi_int], coulomb_spectral
end


"""
    epsilon_causality_vDOSW(electronic_spec, Z_ongrid, phi_ph_ongrid, phi_c, shift_ongrid, w_prime)

Lean vDOS+W ε-integration for the μ-update. Same moment loop as
`epsilon_helpers_vDOSW`, but returns only the Z (causality) and χ (charge)
ω′-integrands. The φ branch and the Coulomb spectral matrix (hence W) are not
needed to conserve charge, so they are dropped — this runs on every root-finder
evaluation of `diff_Ne_realAxis`.
"""
function epsilon_causality_vDOSW(electronic_spec, Z_ongrid, phi_ph_ongrid, phi_c, shift_ongrid, w_prime)
    eps_1, eps_2, dos_1, dos_2, dε, ddos, _, _ = electronic_spec

    M = length(w_prime)
    Nint = length(dε)

    z_int = zeros(Float64, M)
    chi_int = zeros(Float64, M)

    wZ = w_prime .* Z_ongrid

    @inbounds for jε in 1:Nint
        xl = eps_1[jε]
        xr = eps_2[jε]
        invdε = 1.0 / dε[jε]
        m1 = ddos[jε] * invdε
        m0 = dos_1[jε] - xl * m1

        # φ(ω,ε) = φ_ph(ω) + φ_c(ε): the ε-slope Φ1 = φ_c'(ε) is ω-independent
        Φ1 = (phi_c[jε+1] - phi_c[jε]) * invdε
        inv_scale = 1 / (1 + Φ1^2)

        for iw in 1:M
            Φ0 = phi_ph_ongrid[iw] + phi_c[jε] - xl * Φ1
            sh = shift_ongrid[iw]
            wz = wZ[iw]

            S = (sh + Φ0 * Φ1) * inv_scale
            P = sqrt((wz^2 - sh^2 - Φ0^2) * inv_scale + S^2)

            Sp = S + P
            Sm = S - P
            Rp = real(Sp); Ip = imag(Sp)
            Rm = real(Sm); Im_ = imag(Sm)

            i0pp, i1pp, i2pp, i3pp = lorentzian_interval_moments(xl, xr, Rp, Ip)
            i0mm, i1mm, i2mm, i3mm = lorentzian_interval_moments(xl, xr, Rm, Im_)

            c = inv_scale / (2 * P)

            # Z branch: g = ωZ/scale, no slope
            gz = wz * c
            C0, C1 = spectral_C0C1(real(gz), imag(gz), Rp, Rm, Ip, Im_,
                                   i0pp, i1pp, i2pp, i0mm, i1mm, i2mm)
            z_int[iw] += m0 * C0 + m1 * C1

            # χ branch: g = χ/scale, slope = 1/scale
            gchi = sh * c
            C0, C1 = spectral_C0C1(real(gchi), imag(gchi), Rp, Rm, Ip, Im_,
                                   i0pp, i1pp, i2pp, i0mm, i1mm, i2mm)
            D0, D1 = spectral_C0C1(real(c), imag(c), Rp, Rm, Ip, Im_,
                                   i1pp, i2pp, i3pp, i1mm, i2mm, i3mm)
            chi_int[iw] += m0 * (C0 + D0) + m1 * (C1 + D1)
        end
    end

    return z_int, chi_int
end


"""
    lorentzian_moments_weighted(x0, x1, A, B, M0, M1)

Evaluate lorentzian integrals together with weighting factors M0/M1
Allows to sum directly over epsilon
Returns M0_J^(k) in my notation
A, B are R+/I+ or R-/I-
"""
function lorentzian_moments_weighted(x0, x1, A::AbstractVector, B::AbstractVector, M0, M1)
    xl = vec(x0); xr = vec(x1); m0 = vec(M0); m1 = vec(M1)
    Nε = length(xl)
    M  = length(A)
    M0_I0 = zeros(M); M0_I1 = zeros(M); M0_I2 = zeros(M); M0_I3 = zeros(M)
    M1_I0 = zeros(M); M1_I1 = zeros(M); M1_I2 = zeros(M); M1_I3 = zeros(M)
    @inbounds for jε in 1:Nε
        l = xl[jε]; r = xr[jε]; w0 = m0[jε]; w1 = m1[jε]
        @simd for iw in 1:M
            a = A[iw]; b = B[iw]
            b2 = b^2; a2 = a^2
            i0l = atan((a + l) / b) / b;  i0r = atan((a + r) / b) / b
            hl  = log((a + l)^2 + b2);    hr  = log((a + r)^2 + b2)
            i1l = 0.5 * hl - a * i0l;     i1r = 0.5 * hr - a * i0r
            i2l = -a * hl + (a2 - b2) * i0l + l;  i2r = -a * hr + (a2 - b2) * i0r + r
            i3l = i0l * (3 * a * b2 - a * a2) + 0.5 * (3 * a2 - b2) * hl + 0.5 * l * (l - 4 * a)
            i3r = i0r * (3 * a * b2 - a * a2) + 0.5 * (3 * a2 - b2) * hr + 0.5 * r * (r - 4 * a)
            d0 = i0r - i0l; d1 = i1r - i1l; d2 = i2r - i2l; d3 = i3r - i3l
            M0_I0[iw] += w0 * d0; M0_I1[iw] += w0 * d1; M0_I2[iw] += w0 * d2; M0_I3[iw] += w0 * d3
            M1_I0[iw] += w1 * d0; M1_I1[iw] += w1 * d1; M1_I2[iw] += w1 * d2; M1_I3[iw] += w1 * d3
        end
    end

    return (M0_I0, M0_I1, M0_I2, M0_I3, M1_I0, M1_I1, M1_I2, M1_I3)
end


"""
    eval_spectral_integrals(g, M0, M1, ε_p, Rplus, Rminus, Iplus, Iminus, eps_1, eps_2, shift_int; g_slope)

Evaluate spectral integrals for wZ, delta, chi based on the precomputed lorentzian integrals
The lorentzian moments are already weighted by M0/M1 and the ε-sumamtion is done
"""
function eval_spectral_integrals(g, lorentzian_moments_p, lorentzian_moments_m, ε_p, Rplus, Rminus, Iplus, Iminus, shift_int::Bool=false; g_slope=nothing)
    slope = (shift_int && isnothing(g_slope)) ? one.(ε_p) : g_slope

    M0_I0_pp, M0_I1_pp, M0_I2_pp, M0_I3_pp, M1_I0_pp, M1_I1_pp, M1_I2_pp, M1_I3_pp = lorentzian_moments_p
    M0_I0_mm, M0_I1_mm, M0_I2_mm, M0_I3_mm, M1_I0_mm, M1_I1_mm, M1_I2_mm, M1_I3_mm = lorentzian_moments_m

    # spectral coefficients Rg, Ig, R+-,I+-
    C0s, _ = eval_spectral_integral_coefficients(g, ε_p, Rplus, Rminus, Iplus, Iminus,
                                                 M0_I0_pp, M0_I0_mm, M0_I1_pp, M0_I1_mm, M0_I2_pp, M0_I2_mm, M0_I3_pp, M0_I3_mm; g_slope=slope)
    _, C1s = eval_spectral_integral_coefficients(g, ε_p, Rplus, Rminus, Iplus, Iminus,
                                                 M1_I0_pp, M1_I0_mm, M1_I1_pp, M1_I1_mm, M1_I2_pp, M1_I2_mm, M1_I3_pp, M1_I3_mm; g_slope=slope)
    return C0s .+ C1s
end


function eval_spectral_integral_coefficients(g, ε_p, Rplus, Rminus, Iplus, Iminus, I0_pp, I0_mm, I1_pp, I1_mm, I2_pp, I2_mm, I3_pp, I3_mm; g_slope=nothing)
    Rg = real(g./ (2*ε_p))
    Ig = imag(g./ (2*ε_p))

    C0 = @. Ig * (I1_pp - I1_mm)
    C1 = @. Ig * (I2_pp - I2_mm)
    @. C0 += Ig * (Rplus * I0_pp - Rminus * I0_mm)
    @. C1 += Ig * (Rplus * I1_pp - Rminus * I1_mm)
    @. C0 += Rg * (-Iplus * I0_pp + Iminus * I0_mm)
    @. C1 += Rg * (-Iplus * I1_pp + Iminus * I1_mm)

    if !isnothing(g_slope)
        # Add the ε-dependent part of an affine numerator g(ε)=g+g_slope*ε.
        # The vDOS+μ χ integral corresponds to g_slope=1.
        Rg = real(g_slope ./ (2 * ε_p))
        Ig = imag(g_slope ./ (2 * ε_p))

        @. C0 += Ig * (I2_pp - I2_mm)
        @. C1 += Ig * (I3_pp - I3_mm)
        @. C0 += Ig * (Rplus * I1_pp - Rminus * I1_mm)
        @. C1 += Ig * (Rplus * I2_pp - Rminus * I2_mm)
        @. C0 += -Rg * (Iplus * I1_pp - Iminus * I1_mm)
        @. C1 += -Rg * (Iplus * I2_pp - Iminus * I2_mm)
    end

    return C0, C1
end
