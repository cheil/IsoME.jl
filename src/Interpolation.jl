"""
    File containing the interpolation of Weep and Dos onto a uniform grid

    Julia Packages:
        - Interpolations
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


function uniformGrid(start, stop, step)
    if step < 0
        return reverse(start:step:stop)
    else
        return start:step:stop
    end
end




