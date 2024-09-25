"""
    File containing some intergration methods:
        - Gauss Quadrature
        - Simpson rule

    Julia Packages:
        - Interpolations
        - FastGaussQuadrature

    Comments:
        - Gauss-Quadrature maybe not a good choice because:
            a) not guaranteed to include ef in vDOS (energies are not shifted)
            b) roots may lie outside of energy window of Weep/Dos --> extrapolation needed

"""

### Gauss Quadrature ###
function gaussQuad(f, w)
    """
    Perform gauss-quadrature for given weigths and function values at roots
    If f is a matrix perform the integration for each row

    -------------------------------------------------------------------
    Input:
        f:  function values at grid points
        w:  weights 

    --------------------------------------------------------------------
    Output:
        I:  Result  

    --------------------------------------------------------------------
    Comments:
        - 
    --------------------------------------------------------------------
    """


    # 
    if length(f) == length(w)
        I = dot(w, f);
    else
        I = f*w
    end

    return I

end

# jacobi for gauss-Legendre
function jacobi(u, a, b)
    return @. pi / 2 *((tan(pi/2*u))^2 + 1)

    # linear trafo
    # return @. (b-a)/2
end

# inv variable transformation for gauss-Legendre
# transform from u to epsilon
function inv_variable_trafo(u, a, b)
    return @. tan(pi/2*u)

    # linear trafo 
    # return @. ((b-a)*u + a + b)/2
end

### How to use Gauss Quadrature ###
#= 
    # get roots, weights
    roots, weights = gausslegendre(N_Gauss_Quadrature)
    u = roots
    weights .*= jacobi(u)

    epsilon = inv_variable_trafo(u)

    # Calculate DoS at zeros of Legendre Polynoms
    dos_e_data2 = itp_dos(epsilon)  

    # Calculate Weep at zeros of Legendre Polynoms
    Weep2 = itp_Weep(epsilon,epsilon)
=#


### simpson 1/3 rule ### 
function simps(f, a, b)
    """
    Simpson 1/3 integration

    -------------------------------------------------------------------
    Input:
        f:  function values at grid points
        a:  interval start         
        b:  interval end

    --------------------------------------------------------------------
    Output:
        I:  Result  

    --------------------------------------------------------------------
    Comments:
        -
    --------------------------------------------------------------------
    """

    n = length(f) - 1
    h = (b - a) / n
    I = f[1] + f[end]
    f1 = f[2:2:end]
    f2 = f[3:2:end-2]
    for k in eachindex(f1)
        I += 4*f1[k] 
    end
    for k in eachindex(f2)
        I += 2*f2[k] 
    end

    # old version, slightly slower
    # I = h / 3 * (f[1] + 2 * sum(f[3:2:end-2]) + 4 * sum(f[2:2:end]) + f[end])
    
    return I *h/3
end

function simps(x, f)
    """
    Simpson 1/3 integration

    -------------------------------------------------------------------
    Input:
        f:  function values at grid points, MxN matrix
            integration along columns
        x:  grid points, Nx1 vector

    --------------------------------------------------------------------
    Output:
        I:  Result  

    --------------------------------------------------------------------
    Comments:
        -
    --------------------------------------------------------------------
    """


    # integration bounds
    a = x[1]
    b = x[end]

    if length(x) == length(f)
        # stepsize
        n = length(f) - 1
        h = (b - a) / n

        I = h / 3 * (f[1] + 2 * sum(f[3:2:end-2]) + 4 * sum(f[2:2:end]) + f[end])

    else
        # stepsize
        n = size(f,2) - 1
        h = (b - a) / n

        # integrate    
        I = h / 3 * (f[:,1] + 2 * sum(f[:,3:2:end-2], dims=2) + 4 * sum(f[:,2:2:end], dims=2) + f[:,end])
    end


    return I 
end
