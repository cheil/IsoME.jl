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
function bisection(f::Function, a::Number, b::Number, tol::AbstractFloat=1e-6, ftol::AbstractFloat=1e-10, maxiter::Integer=1000)


    # intitial function evaluation
    fa =  f(a)
    fb = f(b)

    fa * fb <= 0 || error("No real root in [a,b]")

    # init
    it = 0
    c = 0
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
function RegulaFalsi(f::Function, a::Number, b::Number, tol::AbstractFloat=1e-10, ftol::AbstractFloat=1e-10,  maxiter::Integer=1000)
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


    # intitial function evaluation
    fa =  f(a)
    fb = f(b)

    fa * fb <= 0 || error("No real root in [a,b]")

    # init
    it = 0
    c = 0
    while abs(b - a) > tol
        # max iterations
        it +=1 
        it != maxiter || error("max number of iterations exceeded")

        # new guess
        c = (a*fb - b*fa) / (fb - fa)

        # function evaluation at new point
        fc = f(c)

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

    return c
end