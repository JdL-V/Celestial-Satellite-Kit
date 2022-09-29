function mars_atm_exp(z)
    T::Float64 = (-23.4 - 0.00222*z)*(z >= 7e3) + (-31 - 0.000998*z)*(z < 7e3);
    T = -170*(T < -170) + T*(T >= -170)
    p::Float64 = 0.699 * exp(-0.00009 * z)
    rho::Float64 = p/(0.1921 * (T + 273.15))
    return rho, T, p
end