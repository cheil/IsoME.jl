"""
    mcmillan_applicable(lambda, mu_star)

Whether the McMillan exponent -1.04(1+λ) / (λ - μ*(1+0.62λ)) is negative. Otherwise
(λ ≲ μ*) it turns positive and the formula returns an arbitrarily large Tc instead of
predicting no superconductivity.
"""
mcmillan_applicable(lambda, mu_star) = lambda - mu_star * (1 + 0.62 * lambda) > 0


"""
    calc_AD_Tc(omega, a2f, mu_star)

Calculate superconducting properties from Allen-Dynes McMillan formula using the energy in meV,
the Eliashberg spectral funciton in 1/meV and the Anderson pseudopotential.

Outside the range of the McMillan formula (see [`mcmillan_applicable`](@ref)) Tc_AD, Tc_ML
and the gap are 0. If the ML correction factor has no finite, positive value (its pole near
λ ≈ 1), Tc_ML falls back to Tc_AD.
"""
function calc_AD_Tc(omega, a2f, mu_star)

    a2f_eps = 1e-2
    a2f   = a2f[omega .> a2f_eps]
    omega = omega[omega .> a2f_eps]

    integrand  = a2f ./ omega 
    lambda     = 2*trapz(omega, integrand)

    integrand  = a2f .*log.(omega) ./ omega 
    omega_log  = 2/lambda * trapz(omega, integrand)
    omega_log  = exp(omega_log)

    # McMillan predicts no superconductivity for λ ≲ μ*: the exponent would turn positive
    Tc = mcmillan_applicable(lambda, mu_star) ?
         omega_log/1.2 * exp( -(1.04*(1+lambda))/(lambda - mu_star - 0.62*lambda*mu_star) ) : 0.0 # Tc in meV

    # omega_1    =   2/lambda*trapz(omega,a2f)
    omega_2    = ( 2/lambda*trapz(omega,a2f.*omega) )^(1/2)

    D1 = 2.46*(1+3.8*mu_star);
    D2 = 1.82*(1+6.3*mu_star)*omega_2./omega_log;

    f1 = ( 1 + (lambda/D1).^(3/2) )^(1/3);
    f2 = 1 + (omega_2/omega_log - 1)*lambda^2 / (lambda^2+D2^2);
    
    Tc_AD = f1*f2*Tc; # Tc in meV
    gap0  = 1.764*Tc_AD

    fom = 1.92 * ( (lambda + omega_log./omega_2 - mu_star.^(1/3) ) ./ (sqrt(lambda) .* exp(omega_log./omega_2) ) ) - 0.08
    fmu = 1.00 + 6.86 .* exp(-lambda./mu_star) ./ (1.0/lambda - mu_star - omega_log./omega_2)
    Tc_ML = fom*fmu*Tc; # Tc in meV
    if !isfinite(fmu) || fmu < 0 || !isfinite(Tc_ML)
        Tc_ML = Tc_AD   # pole of the ML correction factor
    end

    data = [Tc_ML, Tc_AD, gap0, lambda, omega_log]

    return data
end