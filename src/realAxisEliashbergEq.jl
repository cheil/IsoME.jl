"""
file containing the real axis eliashberg equations

@MIT: 
You can add the real axis equations to the functions below 
The functions will be iteratively called in the realAxisSolver.jl file 
Of course if you prefer a different structure for the files/functions just change it as you like
"""


"""

real axis Eliashberg equations in vDOS+W approximation
"""
function realEliashbergEq(itemp::Number, a2f_omega::StepRangeLen{Float64}, a2f::Vector{Float64}, dosef::Float64, ndos::Int64,
                         dos_en::Vector{Float64}, dos::Vector{Float64}, Weep::Matrix{Float64},
                         znormip::Vector{Float64}, phiphip::Vector{Float64}, phicip::Vector{Float64},
                         shiftip::Vector{Float64}, wgCoulomb::Float64)
    """
    vDOS + W equations:
    - itemp: current temperature
    _ a2f_omega: alpha^2F energy grid
    - a2f: alpha^2F values
    - dosef: dos at fermi energy
    - ndos: length of dos vector
    - dos_en: energy grid dos
    - dos: dos values
    - Weep: W(ε,ε')-values, matrix
    - znormip: Z(ω) previous iteration
    - phiphip: ϕ_ph(ω) previous iteration, phonon part
    - phicip: ϕ_c(ω) previous iteration, coulomb part
    - shiftip: χ(ω) previous iterations
    - wcCoulomb: weigth coulomb interaction, damping during first iterations
    """

end


"""

real axis Eliashberg equations in vDOS+μ approximation
"""
function realEliashbergEq()
    # vDOS + mu equations

end



"""

real axis Eliashberg equations in cDOS+μ approximation
"""
function realEliashbergEq(itemp, a2f_omega, a2f, mu, wgCoulomb, deltaip)
    # cDOS + mu equations

end