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
function RegulaFalsi(f::F, fa::Float64, fb::Float64, a::Float64, b::Float64, tol::Float64=1e-6, ftol::Float64=1e-10, maxiter::Int=1000) where {F}
    """
    The Regula Falsi method is a root finding method superior to bisection
    Two function values with opposite sign need to be known

    --------------------------------------------------------------------
    Input:
        f:          Julia function of one argument
        a:          lower boundary
        b:          upper boundary
        tol:        tolerance 
        ftol:       tolerance in f
        maxiter:    maximum number of iterations

    --------------------------------------------------------------------
    Output:
        c:       root

    --------------------------------------------------------------------
    Comments:
        - We are using two different convergence criteria
            * |b-a| < tol
            * |f(root)| < ftol
        
    -------------------------------------------------------------------- 
    """

    fa * fb <= 0 || error("No real root in [a,b]")

    # init
    it = 0
    c = 0.0     # Float64 so the return type is Float64, not Union{Float64,Int64}
    while abs(b - a) > tol
        # max iterations
        it +=1 
        it != maxiter || error("max number of iterations exceeded")

        # new guess
        c = (a*fb - b*fa) / (fb - fa)

        # function evaluation at new point
        fc = f(c)

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

        if abs(fc) < ftol
            break
        elseif fa * fc > 0
            a = c   # Root is in the right half of [a,b].
            fa = fc
        else
            b = c   # Root is in the left half of [a,b].
            fb = fc
        end
    end

    # --------- mu-update progress line ---------
    # erase the progress line (terminal only), cursor back to column 0 so the next
    # table row overwrites it 
    # stdout isa Base.TTY && (print("\r", " "^110, "\r"); flush(stdout))

    return c
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
    phi_shape::Union{Tuple{Int}, NTuple{2,Int}}   # φ is a 1D (cDOS) or 2D (vDOS) array
    started::Bool
end

BroydenMixer(m::Int=8) = BroydenMixer(m, Float64[], Float64[], Vector{Float64}[], Vector{Float64}[], 0, 0, (0, 0), false)

# stack (Z, χ, φ) into one real vector: [Re Z; Im Z; Re χ; Im χ; Re vec(φ); Im vec(φ)]
_broyden_flatten(Z, chi, phi) = vcat(real(Z), imag(Z), real(chi), imag(chi), real(vec(phi)), imag(vec(phi)))

function _broyden_unflatten(mx::BroydenMixer, x::Vector{Float64})
    nz, nc = mx.nz, mx.nc
    np = prod(mx.phi_shape)
    o = 0
    Z   = complex.(x[o+1:o+nz], x[o+nz+1:o+2nz]); o += 2nz
    chi = complex.(x[o+1:o+nc], x[o+nc+1:o+2nc]); o += 2nc
    phi = reshape(complex.(x[o+1:o+np], x[o+np+1:o+2np]), mx.phi_shape)
    return Z, chi, phi
end

"""
    broyden_mix!(mx, β, Z_prev, chi_prev, phi_prev, Z_out, chi_out, phi_out)

One Broyden 2nd-method step. `(*_prev)` is the input to realEliashbergEq, `(*_out)` its
raw output. Returns the mixed `(Z, χ, φ)`. The first call falls back to linear mixing
`x + β(F−x)` and seeds the history.
"""
function broyden_mix!(mx::BroydenMixer, β, Z_prev, chi_prev, phi_prev, Z_out, chi_out, phi_out)
    if !mx.started
        mx.nz = length(Z_prev); mx.nc = length(chi_prev); mx.phi_shape = size(phi_prev)
    end

    xin  = _broyden_flatten(Z_prev, chi_prev, phi_prev)
    xout = _broyden_flatten(Z_out,  chi_out,  phi_out)
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
