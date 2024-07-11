"""
    File containing a linear interpolation of
        - Weep 
        - Dos

    Julia Packages:
        - Interpolations

    Comments:
        - In the new code only the Dos interpolation is of relevance!!

"""

########### Interpolate Weep ###########
function interpolateWeep(epsilon, Weep, interval, intFlag, Nint=nothing)
    """
    !!! THIS FUNCTION IS NOT NEEDED ANYMORE !!!
    
    Interpolate Weep matrix 

    -------------------------------------------------------------------
    Input:
        epsilon:    energy grid
        Weep:       coulomb interaction W(e,e')
        interval:   new energy grid, [interval[1], interval[2]]
                    only relevant for simpson rule
        intFlag:    specify integration method
                        - 0: Gauss-Quadrature
                        - 1: Simpson-Rule
        Nint:      Number of grid points considered in integration scheme
        

    --------------------------------------------------------------------
    Output:
        epsilon:    energies on integration grid
        dos:        dos at epsilon
        Weep:       Coulomb interaction W(e,e')
        roots:      roots used in Gauss-Quadrature 
        weights:    weights used in Gauss-Quadrature
    
    --------------------------------------------------------------------
    Comments:
        - Bug: What if roots are outside of interval when using gauss quad?
    --------------------------------------------------------------------
    """


    ### Interpolation ###
    # Interpolation Object Weep
    epsilonItp = range(epsilon[1], epsilon[end], length(epsilon))       # interpolation vector must be of the form a:b
    itp_Weep = linear_interpolation((epsilonItp, epsilonItp), Weep)
    
    ### new integration grid ###
    if intFlag == 0
        # number of grid points
        if isnothing(Nint)
            # default value gauss
            Nint = 51

        elseif mod(Nint,2) == 0
            # check if Nint is odd to include ef
            # all odd orders of Hermite Polynomials have 0(=ef) as root
            Nint = Nint+1

        end

        # get roots, weights
        roots, weights = gausslegendre(N_Gauss_Quadrature)
        u = roots
        epsilon = inv_variable_trafo(u, epsilonItp[1], epsilonItp[end])


    elseif intFlag == 1
        # number of grid points
        if isnothing(Nint)
            # default value simpson
            Nint = 1000

        # elseif mod(Nint,2) == 0
        #    # simpson is only more exact for an odd number of grid points
        #    Nint = Nint+1

        end

        # define epsilon
        epsilon = collect(range(interval[1], interval[2], Nint))

        # empty roots, weights
        roots = nothing
        weigths = nothing

    else
        error("Invalid integration method! Set intFlag either to 0 or 1!")
    end

    # Calculate Weep at zeros of Legendre Polynoms
    Weep = itp_Weep(epsilon, epsilon)


    return Weep
end


########### Interpolate Dos ###########
function interpolateDos(epsilon, dos, interval, Nint=nothing)
    """
    Read in Weep file to solve the isotropic Migdal-Eliashberg equations

    -------------------------------------------------------------------
    Input:
        epsilon:    energy grid 
        dos:        dos at epsilon
        interval:   new energy grid, [interval[1], interval[2]]
        intFlag:    specify integration method
                        - 0: Gauss-Quadrature
                        - 1: Simpson-Rule
        Nint:      Number of grid points considered in integration scheme
        

    --------------------------------------------------------------------
    Output:
        epsilon:    energies on integration grid
        dos:        dos at epsilon
        roots:      roots used in Gauss-Quadrature 
        weights:    weights used in Gauss-Quadrature
    
    --------------------------------------------------------------------
    Comments:
        -
    --------------------------------------------------------------------
    """


    ### Defaul Values ###
    (~isnothing(Nint)) || (Nint = 1000)
    
    ### Interpolation ###
    # Interpolation Object DoS
    epsilonItp = range(epsilon[1], epsilon[end], length(epsilon))       # interpolation vector must be of the form a:b
    itp_dos = linear_interpolation(epsilonItp, dos)   

    # define epsilon
    epsilon = collect(range(interval[1], interval[2], Nint))

    # Calculate DoS at energy grid points
    dos = itp_dos(epsilon)

    return epsilon, dos
end

