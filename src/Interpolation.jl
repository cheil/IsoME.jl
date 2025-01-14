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
    interpolateWeep(epsilon, Weep,  en_range, bndItp, step)

Linear interpolation of Weep. 
Piecewise interpolation between interval[i] and interval[i+1] 
with grid specifications grid[i].
In grid specify either the step size or the number of grid points.
"""
function interpolateWeep(epsilon, Weep, en_range, bndItp, step)

    ### Interpolation ###
    # Interpolation Object Weep
    epsilonItp = range(epsilon[1], epsilon[end], length(epsilon))       # interpolation vector must be of the form a:b
    itp_Weep = linear_interpolation((epsilonItp, epsilonItp), Weep)
    itp_Weep = scale(interpolate(Weep, BSpline(Constant())), (epsilonItp, epsilonItp)) 
    
    # define new grid
    epsilon = append!(reverse(collect(range(-bndItp-step[2], en_range[1], step = -step[2]))), collect(range(-bndItp, bndItp, step = step[1])),
                    collect(range(bndItp+step[2], en_range[2], step = step[2])))

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
function interpolateDos(epsilon, dos, en_range, bndItp, step)    
    
    ### Interpolation ###
    # Interpolation Object DoS
    epsilonItp = range(epsilon[1], epsilon[end], length(epsilon))       # interpolation vector must be of the form a:b
    itp_dos = scale(interpolate(dos, BSpline(Cubic())), epsilonItp) 
   
    # define new grid
    epsilon = append!(reverse(collect(range(-bndItp-step[2], en_range[1], step = -step[2]))), collect(range(-bndItp, bndItp, step = step[1])),
                    collect(range(bndItp+step[2], en_range[2], step = step[2])))

    # Calculate DoS at energy grid points
    dos = itp_dos(epsilon)

    return epsilon, dos
end


"""
     nonUniformGrid(start_point, end_point; density_factor=1.01, min_step=1, max_step=50)


"""
function nonUniformGrid(start_point, end_point; density_factor=1.01, min_step=1, max_step=50)

# Initialize the grid
grid = [start_point]
current_point = start_point
step = min_step  # Start with the minimum step size

# Generate the grid with increasing spacing
while current_point < end_point
    step = min(step * density_factor, max_step)  # Increment step size but cap it at max_step
    current_point += step
    if current_point <= end_point
        push!(grid, current_point)
    end
end

# Ensure the grid ends exactly at the end_point
if grid[end] < end_point
    push!(grid, end_point)
end

return grid

end

