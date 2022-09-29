function mu = mars_atm_mu(z::Float64)
    T::Float64 = (-23.4 - 0.00222*z)*(z >= 7e3) + (-31 - 0.000998*z)*(z < 7e3);
    T = -170*(T < -170) + T*(T >= -170);
    S::Float64 = 222.
    T0::Float64 = 273.15
    mu0::Float64 = 1.37e-5
    T = T + 273.15;
    mu::Float64 = mu0*(T/T0)^(3/2)*(T0 + S)/(T + S)
    return mu
end