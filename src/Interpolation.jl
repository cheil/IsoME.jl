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
function interpolateWeep(epsilon, Weep, interval, Nint=1000)
    """    
    Interpolate Weep matrix 

    -------------------------------------------------------------------
    Input:
        epsilon:    energy grid
        Weep:       coulomb interaction W(e,e')
        interval:   new energy grid, [interval[1], interval[2]]
                    only relevant for simpson rule
        Nint:      Number of grid points considered in integration scheme
        

    --------------------------------------------------------------------
    Output:
        epsilon:    energies on integration grid
        dos:        dos at epsilon
        Weep:       Coulomb interaction W(e,e')
    
    --------------------------------------------------------------------
    Comments:
        - 
    --------------------------------------------------------------------
    """


    ### Interpolation ###
    # Interpolation Object Weep
    epsilonItp = range(epsilon[1], epsilon[end], length(epsilon))       # interpolation vector must be of the form a:b
    itp_Weep = linear_interpolation((epsilonItp, epsilonItp), Weep)
    itp_Weep = scale(interpolate(Weep, BSpline(Constant())), (epsilonItp, epsilonItp)) 
    
    # define epsilon
    epsilon = collect(range(interval[1], interval[2], Nint))


    # Calculate Weep at zeros of Legendre Polynoms
    Weep = itp_Weep(epsilon, epsilon)


    return Weep
end


########### Interpolate Dos ###########
function interpolateDos(epsilon, dos, interval, Nint=1000)
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
    
    ### Interpolation ###
    # Interpolation Object DoS
    epsilonItp = range(epsilon[1], epsilon[end], length(epsilon))       # interpolation vector must be of the form a:b
    #itp_dos = linear_interpolation(epsilonItp, dos)  
    itp_dos = scale(interpolate(dos, BSpline(Cubic())), epsilonItp) 

    # define epsilon
    epsilon = collect(range(interval[1], interval[2], Nint))

    # Calculate DoS at energy grid points
    dos = itp_dos(epsilon)

    return epsilon, dos
end

