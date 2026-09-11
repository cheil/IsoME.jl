"""
    File containing some mixing (root finding) schemes
        - Broyden
            * function can be a blackbox, i.e. needs just function values
        - Bisection
            * function must be known

    In Self-consistent calculations a mixing scheme is ultimately a root finding routine
    as we are searching for solutions to the problem f(x) = x --> F(x) = f(x) - x = 0

    Julia Packages:
        - 

    Comments:
        - There is also the Julia-Package "Roots" which contains some root finding methods
          but it is slower than our methods

"""


### Bisection
"""
    bisection(f, a, b, tol=1e-6, ftol=1e-10, maxiter=1000)

The bisection method is a simple root-finding method
Two function values with opposite sign need to be known
"""
function bisection(f::F, a::Number, b::Number, tol::Float64=1e-6, ftol::Float64=1e-10, maxiter::Int=1000) where {F}


    # intitial function evaluation
    fa =  f(a)
    fb = f(b)

    fa * fb <= 0 || error("No real root in [a,b]")

    # init
    it = 0
    c = 0.0     # Float64 so the return type is Float64, not Union{Float64,Int64}
    while abs(b - a) > tol
        # max iterations
        it += 1
        it != maxiter || error("max number of iterations exceeded")

        # new guess
        c = (a + b) / 2

        # function evaluation at midpoint
        fc = f(c)

        if abs(fc) < ftol
            break
        elseif fa * fc > 0
            a = c   # Root is in the right half of [a,b].
            fa = fc
        else
            b = c   # Root is in the left half of [a,b].
        end
    end

    return c
end



### Regula Falsi
"""
The Regula Falsi method is a root finding method superior to bisection
Two function values with opposite sign need to be known

Illinois damping is applied only once the same endpoint has been carried over
twice in a row. Plain regula falsi can otherwise leave one endpoint stuck forever
on a convex/concave function, so |b-a| never reaches tol and the iteration count
explodes (exp(x)-5 on [0,20] does not converge in 1000 iterations without the
damping, versus ~20 with it). Damping unconditionally would instead spoil the
well-behaved case: on a smooth monotone f the undamped secant step lands almost
exactly on the root, and halving perturbs it, so the ftol exit triggers at a
noticeably worse point. Waiting for a repeated retention keeps the fast path
untouched and only intervenes when stagnation actually sets in.

"""
function RegulaFalsi(f::F, fa::Float64, fb::Float64, a::Float64, b::Float64, tol::Float64=1e-6, ftol::Float64=1e-10, maxiter::Int=1000) where {F}


    fa * fb <= 0 || error("No real root in [a,b]")

    # init
    it = 0
    # Float64 so the return type is Float64, not Union{Float64,Int64}. Seeded with the
    # better endpoint rather than 0.0, so a bracket that is already narrower than tol
    # returns a point inside it instead of the origin.
    c = abs(fa) <= abs(fb) ? float(a) : float(b)
    kept = 0        # which endpoint was carried over last: 1 = a, 2 = b
    keptRun = 0     # how many iterations in a row it has been carried over
    while abs(b - a) > tol
        # max iterations
        it +=1
        it != maxiter || error("max number of iterations exceeded")

        # new guess; fall back to the midpoint if the secant step is degenerate
        # (fb == fa) or lands outside the bracket, which plain regula falsi can do
        c = fb == fa ? (a + b) / 2 : (a*fb - b*fa) / (fb - fa)
        (min(a, b) < c < max(a, b)) || (c = (a + b) / 2)

        # function evaluation at new point
        fc = f(c)

        # a non-finite value carries no sign information: bisect instead
        if !isfinite(fc)
            c = (a + b) / 2
            fc = f(c)
            isfinite(fc) || error("function is not finite inside the bracket")
        end

        # --------- mu-update progress line --------- 
        # in-place progress. Only on a real terminal (Base.TTY): 
        # there "\r"+erase blanks
        # the line and the following table row (stdout, same cursor) overwrites it cleanly.
        # In Jupyter/IJulia there is no ANSI cursor control and committed lines can't be
        # reclaimed, so we stay silent to avoid corrupting / shifting the table.
        # if stdout isa Base.TTY
        #     msg = "  μ-update: it=$it  |b-a|=$(round(abs(b - a), sigdigits=3)) (tol=$tol)  |fc|=$(round(abs(fc), sigdigits=3)) (ftol=$ftol)"
        #     print("\r", rpad(msg, 110)); flush(stdout)
        # end

        if abs(fc) <= ftol
            break
        elseif fa * fc > 0
            a = c   # Root is in the right half of [a,b].
            fa = fc
            keptRun = kept == 2 ? keptRun + 1 : 1
            kept = 2
            keptRun >= 2 && (fb /= 2)   # Illinois damping, see the docstring
        else
            b = c   # Root is in the left half of [a,b].
            fb = fc
            keptRun = kept == 1 ? keptRun + 1 : 1
            kept = 1
            keptRun >= 2 && (fa /= 2)
        end
    end

    # --------- mu-update progress line ---------
    # erase the progress line (terminal only), cursor back to column 0 so the next
    # table row overwrites it 
    # stdout isa Base.TTY && (print("\r", " "^110, "\r"); flush(stdout))

    return c
end


# ----------------------------------------------------------------
#                   Minimal scalar root finder.
# ----------------------------------------------------------------
# Drop-in replacement for the one Roots.jl entry point IsoME uses:
#   find_zero(f, x0)
#
# Roots' default for a scalar starting point is Order0(), a secant/bisection
# hybrid. The same two-stage strategy is used here: grow a bracket outwards from
# x0, then hand it to RegulaFalsi in Mixing.jl - the bracketed solve is not
# duplicated here. As with Roots, a failure to converge throws; every call site is
# already wrapped in a try/catch with its own fallback, so that behaviour is kept.
#


"""
    find_zero(f, x0; xtol, maxIter) -> Float64

Find a root of the scalar function `f` starting from `x0`.

--------------------------------------------------------------------
Input:
    f:      scalar function of one variable
    x0:     starting point (not a bracket)

Keywords:
    xtol:       relative width of the bracket at which to stop
    maxIter:    maximum iterations of the bracketed solve
    maxExpand:  maximum bracket-expansion steps

--------------------------------------------------------------------
Output:
    the root, as a Float64

--------------------------------------------------------------------
Comments:
    - Throws if no sign change can be found around `x0`, mirroring Roots'
      behaviour on a failed solve, so the callers' fallbacks still trigger.
    - `f` may be undefined away from the root (the fitted models used here
      contain `log`/`sqrt` of quantities that go negative). Points where `f`
      throws or is non-finite are skipped during the expansion instead of
      aborting the search.
"""
function find_zero(f::F, x0::Real; xtol::Float64 = 1e-12, maxIter::Int64 = 200,
                   maxExpand::Int64 = 100) where {F}

    x = float(x0)
    fx = evalSafe(f, x)
    isfinite(fx) || error("find_zero: f is not finite at the starting point x0 = $x0")
    fx == 0 && return x

    # grow a bracket outwards from x0, alternating sides
    step = max(abs(x), 1.0) * 0.05
    for _ in 1:maxExpand
        for s in (1.0, -1.0)
            c = x + s * step
            fc = evalSafe(f, c)
            if isfinite(fc) && fx * fc <= 0
                # hand the bracket to the shared solver. f is wrapped so a model that
                # is undefined somewhere inside the bracket yields NaN rather than
                # throwing, which RegulaFalsi answers with a bisection step.
                tol = xtol * max(1.0, abs(x), abs(c))
                return RegulaFalsi(y -> evalSafe(f, y), fx, fc, x, c, tol, 0.0, maxIter)
            end
        end
        step *= 1.6
        step > 1e13 && break
    end

    error("find_zero: no sign change found around x0 = $x0")
end


"""
    evalSafe(f, x)

Evaluate `f(x)`, mapping a domain error onto `NaN` so the caller can skip the
point instead of failing.
"""
function evalSafe(f::F, x::Float64) where {F}
    try
        v = f(x)
        return isa(v, Real) ? float(v) : NaN
    catch ex
        ex isa InterruptException && rethrow(ex)
        return NaN
    end
end


############################################################
# -------------------- Broyden mixing -------------------- #
############################################################
"""
    BroydenMixer(m)

State for the limited-memory Broyden (2nd / "bad" method) self-consistency mixer.
The fields (Z, χ, φ) are flattened into one real vector (real/imag stacked) so the
quasi-Newton update is done on a single unambiguous residual g(x) = F(x) - x.
`m` is the history depth. β (the underlying linear-mixing weight, = initial inverse
Jacobian −βI) is supplied per call so the Coulomb ramp-up still applies.
"""
mutable struct BroydenMixer
    m::Int
    xin_prev::Vector{Float64}
    F_prev::Vector{Float64}
    u::Vector{Vector{Float64}}      # stored update vectors:  Gₙ = −βI + Σ uᵢ vᵢᵀ
    v::Vector{Vector{Float64}}
    nz::Int                         # component sizes (for un-flatten)
    nc::Int
    n_ph::Int                       # length of φ_ph(ω)
    n_c::Int                        # length of φ_c(ε); 0 in vDOS+μ
    started::Bool
end

BroydenMixer(m::Int=8) = BroydenMixer(m, Float64[], Float64[], Vector{Float64}[], Vector{Float64}[], 0, 0, 0, 0, false)

# stack (Z, χ, φ_ph, φ_c) into one real vector: [Re Z; Im Z; Re χ; Im χ; Re φ_ph; Im φ_ph; Re φ_c; Im φ_c]
_broyden_flatten(Z, chi, phi_ph, phi_c) = vcat(real(Z), imag(Z), real(chi), imag(chi),
                                               real(phi_ph), imag(phi_ph), real(phi_c), imag(phi_c))

function _broyden_unflatten(mx::BroydenMixer, x::Vector{Float64})
    nz, nc, nph, ncphi = mx.nz, mx.nc, mx.n_ph, mx.n_c
    o = 0
    Z      = complex.(x[o+1:o+nz], x[o+nz+1:o+2nz]);      o += 2nz
    chi    = complex.(x[o+1:o+nc], x[o+nc+1:o+2nc]);      o += 2nc
    phi_ph = complex.(x[o+1:o+nph], x[o+nph+1:o+2nph]);   o += 2nph
    phi_c  = complex.(x[o+1:o+ncphi], x[o+ncphi+1:o+2ncphi])
    return Z, chi, phi_ph, phi_c
end

"""
    broyden_mix!(mx, β, Z_prev, chi_prev, phi_ph_prev, phi_c_prev, Z_out, chi_out, phi_ph_out, phi_c_out)

One Broyden 2nd-method step. `(*_prev)` is the input to realEliashbergEq, `(*_out)` its
raw output. Returns the mixed `(Z, χ, φ_ph, φ_c)`. The first call falls back to linear mixing
`x + β(F−x)` and seeds the history.
"""
function broyden_mix!(mx::BroydenMixer, β, Z_prev, chi_prev, phi_ph_prev, phi_c_prev, Z_out, chi_out, phi_ph_out, phi_c_out)
    if !mx.started
        mx.nz = length(Z_prev); mx.nc = length(chi_prev); mx.n_ph = length(phi_ph_prev); mx.n_c = length(phi_c_prev)
    end

    xin  = _broyden_flatten(Z_prev, chi_prev, phi_ph_prev, phi_c_prev)
    xout = _broyden_flatten(Z_out,  chi_out,  phi_ph_out,  phi_c_out)
    F    = xout .- xin                       # residual g(x) = F(x) - x

    if !mx.started
        xnew = xin .+ β .* F                 # first step: plain linear mixing
        mx.xin_prev = xin; mx.F_prev = F; mx.started = true
        return _broyden_unflatten(mx, xnew)
    end

    # rank-1 update of the inverse Jacobian (Broyden II), limited memory
    dX   = xin .- mx.xin_prev
    dF   = F   .- mx.F_prev
    dFdF = dot(dF, dF)
    if dFdF > 1e-30                          # skip update if the residual barely moved (converged)
        un = dX .+ β .* dF
        for i in eachindex(mx.u)
            un .-= mx.u[i] .* dot(mx.v[i], dF)
        end
        push!(mx.u, un)
        push!(mx.v, dF ./ dFdF)
        if length(mx.u) > mx.m               # roll the history window
            popfirst!(mx.u); popfirst!(mx.v)
        end
    end

    # quasi-Newton step:  x_new = xin + β F − Σ uᵢ (vᵢ·F)
    xnew = xin .+ β .* F
    for i in eachindex(mx.u)
        xnew .-= mx.u[i] .* dot(mx.v[i], F)
    end

    mx.xin_prev = xin; mx.F_prev = F
    return _broyden_unflatten(mx, xnew)
end
