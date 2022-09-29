function mu = earth_atm_mu(z)
    if z <= 11e3
        Ti::Float64 = 288.15
        lambda::Float64 = -6.5e-3
        zi::Float64 = 0.
    elseif z <= 20e3
        Ti::Float64 = 216.65
        lambda::Float64 = 0
        zi::Float64 = 11e3
    elseif z <= 32e3
        Ti::Float64 = 216.65
        lambda::Float64 = 1e-3
        zi::Float64 = 20e3
    elseif z <= 47e3
        Ti::Float64 = 228.65
        lambda::Float64 = 2.8e-3
        zi::Float64 = 32e3
    elseif z <= 51e3
        Ti::Float64 = 270.65
        lambda::Float64 = 0
        zi::Float64 = 47e3
    elseif z <= 71e3
        Ti::Float64 = 270.65
        lambda::Float64 = -2.8e-3
        zi::Float64 = 51e3
    else
        Ti::Float64 = 214.65
        lambda::Float64 = -2e-3
        zi::Float64 = 71e3
    end
    T::Float64 = Ti + lambda*(z - zi)
    S::Float64 = 120.
    T0::Float64 = 219.15
    mu0::Float64 = 1.37*10^(-5)
    mu = mu0*(T/T0)^(3/2)*(T0 + S)/(T + S)
    return mu
end
