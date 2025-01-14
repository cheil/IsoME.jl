"""
    File containing a linear interpolation of
        - Weep 
        - Dos

    Julia Packages:
        - Interpolations

    Comments:
        - In the new code only the Dos interpolation is of relevance!!

"""


"""
    interpolateDos(epsilon, dos, interval, grid=[("", 1000)])

Linear interpolation of Weep. 
Piecewise interpolation between interval[i] and interval[i+1] 
with grid specifications grid[i].
In grid specify either the step size or the number of grid points.
"""
function interpolateWeep(epsilon, Weep, interval, grid=[("", 1000)])

    ### check input
    if length(interval)-1 != length(grid)
        error("Number of intervals and grid specifications do not match!")
    end

    ### Interpolation ###
    # Interpolation Object Weep
    epsilonItp = range(epsilon[1], epsilon[end], length(epsilon))       # interpolation vector must be of the form a:b
    itp_Weep = linear_interpolation((epsilonItp, epsilonItp), Weep)
    itp_Weep = scale(interpolate(Weep, BSpline(Constant())), (epsilonItp, epsilonItp)) 
    
    # define new grid
    epsilon = Vector{Float64}()
    for k in eachindex(grid)
        if isempty(grid[k][1]) || grid[k][1] == "N" || grid[k][1] == "length"
            # define epsilon
            append!(epsilon, collect(range(interval[k], interval[k+1], length=grid[k][2])))
        elseif grid[k][1] == "de" || grid[k][1] == "step"
            append!(epsilon, collect(range(interval[k], interval[k+1], step = grid[k][2])))
        end
    end


    # Calculate Weep at zeros of Legendre Polynoms
    Weep = itp_Weep(epsilon, epsilon)


    return Weep
end


"""
    interpolateDos(epsilon, dos, interval, grid=[("", 1000)])

Linear interpolation of dos and energies. 
Piecewise interpolation between interval[i] and interval[i+1] 
with grid specifications grid[i].
In grid specify either the step size or the number of grid points.
"""
function interpolateDos(epsilon, dos, interval, grid=[("", 1000)])

    ### check input
    if length(interval) - 1 != length(grid)
        error("Number of intervals and grid specifications do not match!")
    end
    
    
    ### Interpolation ###
    # Interpolation Object DoS
    epsilonItp = range(epsilon[1], epsilon[end], length(epsilon))       # interpolation vector must be of the form a:b
    #itp_dos = linear_interpolation(epsilonItp, dos)  
    itp_dos = scale(interpolate(dos, BSpline(Cubic())), epsilonItp) 
    # define new grid

   
    epsilon = Vector{Float64}()
    for k in eachindex(grid)
        if isempty(grid[k][1]) || grid[k][1] == "N" || grid[k][1] == "length"
            # define epsilon
            append!(epsilon, collect(range(interval[k], interval[k+1], length=grid[k][2])))
        elseif grid[k][1] == "de" || grid[k][1] == "step"
            append!(epsilon, collect(range(interval[k], interval[k+1], step = grid[k][2])))
        end
    end


    # Calculate DoS at energy grid points
    dos = itp_dos(epsilon)

    return epsilon, dos
end

