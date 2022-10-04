function mu = earth_atm_mu(z)
if z <= 11e3
    Ti = 288.15;
    lambda = -6.5e-3;
    zi = 0;
elseif z <= 20e3
    Ti = 216.65;
    lambda = 0;
    zi = 11e3;
elseif z <= 32e3
    Ti = 216.65;
    lambda = 1e-3;
    zi = 20e3;
elseif z <= 47e3
    Ti = 228.65;
    lambda = 2.8e-3;
    zi = 32e3;
elseif z <= 51e3
    Ti = 270.65;
    lambda = 0;
    zi = 47e3;
elseif z <= 71e3
    Ti = 270.65;
    lambda = -2.8e-3;
    zi = 51e3;
else
    Ti = 214.65;
    lambda = -2e-3;
    zi = 71e3;
end
    T = Ti + lambda*(z-zi);
    S = 120;
    T0 = 219.15;
    mu0 = 1.37*10^(-5);
    mu = mu0*(T/T0)^(3/2)*(T0+S)/(T+S);
end