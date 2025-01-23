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
    interpolateInputs(itpDos, dos_en, itpStepSize, itpBounds, encut; itpWeep=nothing, Wen=nothing)

Linear interpolation of inputs (energies, dos, Weep). 
Piecewise interpolation for intervals in itpBounds with steps itpStepSize
"""
function interpolateInputs(itpDos, dos_en, itpStepSize, itpBounds, encut; itpWeep=nothing, Wen=nothing)

    if (length(itpStepSize)-1) != length(itpBounds)
        error("Number of interpolation steps and interpolation bounds do not match!")
    end

    # set up interval
    bndItp = abs.(itpBounds)
    steps = append!(-reverse(itpStepSize), itpStepSize)
    if isnothing(Wen)
        enStart = [dos_en[1]]
        enEnd = [dos_en[end]]
    else
        enStart = [Wen[findfirst(Wen .> dos_en[1])]]
        enEnd = [Wen[findlast(Wen .< dos_en[end])]]
    end
    if encut != -1
        enStart = [maximum(push!(enStart, -encut))]
        enEnd = [minimum(push!(enEnd, encut))]
    end

    en_range = append!(enStart, -reverse(bndItp), [0,0], bndItp, enEnd)
    lenRange = length(en_range)
    numIntervals = length(bndItp)
    # check if increasing order
    [(en_range[i] > en_range[i-1]) || (en_range[i] = en_range[i-1]) for i in 2:(numIntervals+1)]
    [(en_range[i] < en_range[i+1]) || (en_range[i] = en_range[i+1]) for i in lenRange-1:-1:(lenRange-numIntervals)]

    # (start, stop, step)
    interval = tuple.(en_range[2:end-1], en_range[vcat(1:(1+numIntervals), ((lenRange-numIntervals):lenRange))], steps)

    # define new grid
    epsilon = Vector{Float64}()
    for (start, stop, step) in interval
        append!(epsilon, uniformGrid(start, stop, step))
        #append!(epsilon, nonUniformGrid(start, stop; density_factor=1.05, min_step=step, max_step=50*step))
     end
     epsilon = unique(epsilon)
 
     # Calculate DoS at energy grid points
     dos = itpDos(epsilon)

     if ~isnothing(itpWeep)
        # Calculate Weep at energy grid points
        Weep = itpWeep(epsilon, epsilon)
     else
        Weep = nothing
     end
 
     return epsilon, dos, Weep

end

"""
    interpolateWeep(epsilon, Weep,  interval)

Linear interpolation of Weep. 
Piecewise interpolation for interval start:step:stop
"""
function interpolateWeep(epsilon, Weep, interval)

    ### Interpolation ###
    # Interpolation Object Weep
    epsilonItp = range(epsilon[1], epsilon[end], length(epsilon))       # interpolation vector must be of the form a:b
    itp_Weep = scale(interpolate(Weep, BSpline(Constant())), (epsilonItp, epsilonItp))

    # define new grid
    epsilon = Vector{Float64}()
    for (start, stop, step) in interval
        append!(epsilon, uniformGrid(start, stop, step))
    end
    epsilon = unique(epsilon)

    # Calculate Weep at zeros of Legendre Polynoms
    Weep = itp_Weep(epsilon, epsilon)

    return Weep
end


"""
    interpolateDos(epsilon, dos, interval)

Linear interpolation of dos and energies. 
Piecewise interpolation between interval[i] and interval[i+1] 
with grid specifications grid[i].
In grid specify either the step size or the number of grid points.
"""
function interpolateDos(epsilon, dos, interval)  

    ### Interpolation ###
    # Interpolation Object DoS
    epsilonItp = range(epsilon[1], epsilon[end], length(epsilon))       # interpolation vector must be of the form a:b
    itp_dos = scale(interpolate(dos, BSpline(Cubic())), epsilonItp)

    # define new grid
    epsilon = Vector{Float64}()
    for (start, stop, step) in interval
        append!(epsilon, uniformGrid(start, stop, step))
        #append!(epsilon, nonUniformGrid(start, stop; density_factor=1.05, min_step=step, max_step=50*step))
    end
    epsilon = unique(epsilon)

    # Calculate DoS at energy grid points
    dos = itp_dos(epsilon)

    return epsilon, dos
end


function uniformGrid(start, stop, step)
    if step < 0
        return reverse(start:step:stop)
    else
        return start:step:stop
    end
end


"""
     nonUniformGrid(start_point, end_point; density_factor=1.01, min_step=1, max_step=50)


"""
function nonUniformGrid(start, stop; density_factor=1.05, min_step=1, max_step=50)

# Initialize the grid
grid = Vector{Float64}([start])
current_point = start
step = min_step  # Start with the minimum step size

# Generate the grid with increasing spacing
while abs(current_point) < abs(stop)
    if step <0
        step = max(step * density_factor, max_step)  # Increment step size but cap it at max_step
    else
        step = min(step * density_factor, max_step)  # Increment step size but cap it at max_step
    end
    current_point += step
    if abs(current_point) <= abs(stop)
        push!(grid, current_point)
    end
end

# Ensure the grid ends exactly at the end_point
if grid[end] < stop
    push!(grid, stop)
end

if step < 0
    return reverse(round.(grid))
else
    return round.(grid)
end

end

