#
# Minimal Levenberg-Marquardt least-squares fit.
#
# Drop-in replacement for the subset of LsqFit.jl that IsoME uses:
#   fit = curve_fit(model, xdata, ydata, p0); par = fit.param
#
# LsqFit's curve_fit defaults to autodiff=:finite, and none of the IsoME call
# sites override it, so the Jacobian is built by forward differences here too -
# the numerics of the original code path are reproduced, not approximated.
#


"""
    FitResult

Result of [`curve_fit`](@ref). `param` holds the fitted parameters, `resid` the
final residual vector `model(x, param) - y`.
"""
struct FitResult
    param::Vector{Float64}
    resid::Vector{Float64}
    converged::Bool
    iterations::Int64
end


"""
    curve_fit(model, xdata, ydata, p0; kwargs...) -> FitResult

Fit `model(xdata, p)` to `ydata` in the least-squares sense, starting from `p0`,
using Levenberg-Marquardt with a forward-difference Jacobian.

--------------------------------------------------------------------
Input:
    model:      callable model(x, p) returning a vector like ydata
    xdata:      independent variable
    ydata:      observations
    p0:         initial parameter guess

Keywords:
    maxIter:    maximum LM iterations
    ftol:       stop when the relative decrease of the residual sum drops below
    xtol:       stop when the relative parameter step drops below
    lambda0:    initial damping
    lambdaUp:   damping growth factor on a rejected step
    lambdaDown: damping shrink factor on an accepted step

--------------------------------------------------------------------
Output:
    FitResult with fields param, resid, converged, iterations

--------------------------------------------------------------------
Comments:
    - The models fitted in IsoME are only defined on part of the parameter
      space (`log`, `sqrt` of quantities that go negative), so a trial step may
      throw or return non-finite values. Such a step is treated as a rejection
      and damping is increased, which keeps the search inside the valid region
      instead of propagating a NaN.
    - The Tc-search fit has 2 data points and 3 parameters, i.e. it is
      underdetermined and J'J is singular. Marquardt's diagonal scaling plus a
      small absolute floor keeps the normal equations solvable; any of the
      infinitely many exact solutions is acceptable there, since the result is
      only used to pick the next (integer) temperature.
"""
function curve_fit(model::F, xdata, ydata, p0::AbstractVector;
                   maxIter::Int64 = 1000, ftol::Float64 = 1e-10,
                   xtol::Float64 = 1e-10, lambda0::Float64 = 1e-3,
                   lambdaUp::Float64 = 10.0, lambdaDown::Float64 = 10.0) where {F}

    p = convert(Vector{Float64}, collect(p0))
    n = length(p)

    r = residual(model, xdata, ydata, p)
    all(isfinite, r) || error("curve_fit: model is not finite at the initial guess p0 = $p0")
    S = sum(abs2, r)

    m = length(r)
    J = zeros(Float64, m, n)
    lambda = lambda0
    converged = false
    it = 0

    while it < maxIter
        it += 1

        # forward-difference Jacobian, same scheme as FiniteDiff's :forward
        for j in 1:n
            h = sqrt(eps(Float64)) * max(abs(p[j]), 1.0)
            pj = copy(p)
            pj[j] += h
            h = pj[j] - p[j]                    # exact representable step
            rj = residual(model, xdata, ydata, pj)
            all(isfinite, rj) || (rj = r)       # flat column rather than NaN
            @. J[:, j] = (rj - r) / h
        end

        JtJ = J' * J
        Jtr = J' * r
        diagScale = max.(diag(JtJ), 1e-12)      # Marquardt scaling with a floor

        # inner loop: grow damping until the step is accepted
        stepTaken = false
        for _ in 1:60
            delta = try
                -((JtJ + lambda * Diagonal(diagScale)) \ Jtr)
            catch ex
                ex isa InterruptException && rethrow(ex)
                fill(NaN, n)
            end

            if all(isfinite, delta)
                pNew = p .+ delta
                rNew = residual(model, xdata, ydata, pNew)

                if all(isfinite, rNew)
                    SNew = sum(abs2, rNew)
                    if SNew < S
                        # accept
                        relS = (S - SNew) / max(S, eps(Float64))
                        relP = norm(delta) / max(norm(p), xtol)
                        p, r, S = pNew, rNew, SNew
                        lambda = max(lambda / lambdaDown, 1e-12)
                        stepTaken = true
                        converged = relS < ftol || relP < xtol
                        break
                    end
                end
            end

            lambda *= lambdaUp
            lambda > 1e14 && break              # damping exhausted
        end

        if !stepTaken
            # no damping value gave a decrease: we sit at a stationary point of
            # the residual sum, which is convergence, not failure
            converged = true
            break
        end
        converged && break
    end

    return FitResult(p, r, converged, it)
end


"""
    residual(model, xdata, ydata, p)

Evaluate `model(xdata, p) - ydata`, mapping a domain error (`log`/`sqrt` of a
negative argument, etc.) onto a non-finite residual so the caller can reject the
step instead of failing.
"""
function residual(model::F, xdata, ydata, p::Vector{Float64}) where {F}
    try
        return collect(Float64, model(xdata, p) .- ydata)
    catch ex
        ex isa InterruptException && rethrow(ex)
        return fill(NaN, length(ydata))
    end
end
