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


### Broyden mixing
function broyden_mixing(iter, x, dx, F, dF, G=I, G0=I, wg0=1e-2, wg=1)
    """
    Modified Broydens method for accelerating self-consistent 
    calculations as suggested in Vanderbilt and Louie, Phys.Rev.B 30, 6118 (1984)
    Relevant equations: (A6), (A14), (A15), (A16)

    Old version:
    calculations as suggested in D.D.Johnson, PRB 38, 12807 (1988)
    Relevant equations: (13a), (13b), (15a), (15b)

    m: number of iterations used in the broyden scheme
    M: length of mixing variable x

    -------------------------------------------------------------------
    Input:
        x:      current guess                               Mx1 matrix
        dx:     x(n+1) - x(n) of the m previous iterations  mxM matrix
        F:      F(x) = 0 of current iteration               Mx1 matrix
                For fixed point search pass F(x) = f(x) - x 
                to this function
        dF:     F(n+1) - F(n) of the m previous iterations  mxM matrix
        u:      G0*dF' + dx' of the m previous iterations   Mxm matrix
        wg0:    weight inverse jacobian                     scalar
        wg:     weight of each previous iteration           mx1 matrix
        G:      updated guess for inverse of jacobian       MxM matrix
        G0:     initial guess for inverse of jacobian       MxM matrix

    --------------------------------------------------------------------
    Output:
        x_new:  updated guess 
        dx:     updated difference of x
        dF:     updated difference of F
        wg:     updated weigths
        G:      updated jacobian
    
    --------------------------------------------------------------------
    Comments:
    When calling this function you need to specify
        - Initial guess x_0 and G0 (suggestion: G0 = I)
        - Initial weigths wg and wg0
        - Initialize dx = zeros(m,M), dx[end,:] = x_0 
            --> dx .= 0 in first iteration
        - Initialize dF = zeros(m,M)          
            --> dF[end, :] = F(x_0) = f(x_0) - x_0 in first iteration
        - The size of dx/dF along dimension 1 defines how many iterations 
          are considered in the Broyden scheme

    --------------------------------------------------------------------
    """

    
    ### Update Input ### 
    if iter > 1
        # dF = F(current) - F(previous)
        dF[1, :] = F - dF[1, :] 

        # dx = x(current) - x(previous)
        dx[1, :] = x - dx[1,:]

        # normalize
        inv_norm = inv(norm(dF[1,:]))
        dF[1,:] = dF[1,:]*inv_norm
        dx[1, :] = dx[1,:]*inv_norm
    end

    #=
    ### Broyden update ###
    # a , mxm matrix
    a = Symmetric(wg * transpose(wg) .* transpose(dF * transpose(dF)))

    # beta , mxm matrix
    beta_broyden = inv(Symmetric(wg0^2 * I + a))

    # u , Mxm matrix
    u = G0 * transpose(dF) + transpose(dx)

    # gamma , mx1 vector
    c_km = wg .* (dF * F)
    gamma = beta_broyden * c_km

    # broyden update
    x_new =  x + G0 * F - u * (wg .* gamma)
    =#


    ### Even another broyden ### 
    # Vanderbilt and Louie, Phys.Rev.B 30, 6118 (1984)
    m = size(dx, 1)
    M = size(dx,2)
    sigma = zeros(m,1)
    beta_broyden = zeros(M,M)
    for k = 1:m
        sigma = 1 + transpose(dx[k,:])*dx[k,:]/wg0^2
        beta_broyden = beta_broyden  - dx[k,:]*transpose(dx[k,:])/sigma/wg0^4
    end
    beta_broyden = beta_broyden + I/(wg0^2)

    gamma = wg0^2*G0 - transpose(wg.*dF)*dx

    if iter == 1
        x_new = x + G0*F
    else
        x_new = x + G*F
    end
    G = inv(gamma*beta_broyden)


    ### Input for next iteration ###
    # write F(previous) to dF for next iteration
    dF[2:end, :] = dF[1:end-1, :]
    dF[1, :] = F

    # write x(previous) to dx for next iteration
    dx[2:end, :] = dx[1:end-1, :]
    dx[1, :] = x

    # update weights, not recommended
    #wg[2:end, :]= wg[1:end-1, :]
    #wg[1] = inv(sqrt(sum(F.*F)))


    # convert 1x1 matrix to scalar
    if length(x_new) == 1
        x_new = x_new[1]
    end
    


    #return x_new, dx, dF, wg     # G not used in first broyden method

    # other broyden return G/G0
    return x_new, dx, dF, G, wg



end


### Bisection
function bisection(f::Function, a::Number, b::Number, tol::AbstractFloat=1e-6, ftol::AbstractFloat=1e-10, maxiter::Integer=1000)
    """
    The bisection method is a simple root-finding method
    Two function values with opposite sign need to be known

    --------------------------------------------------------------------
    Input:
        f:          Julia function of one argument
        a:          lower boundary
        b:          upper boundary
        tol:        tolerance 
        maxiter:    maximum number of iterations

    --------------------------------------------------------------------
    Output:
        c:       root

    --------------------------------------------------------------------
    Comments:
        -
        
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