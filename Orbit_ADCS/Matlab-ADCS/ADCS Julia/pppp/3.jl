# kinetic energy from photons
# ϕ*A/mc*(1 + q)
ϕ = 1e3
A = 1
m = 10
c = 299792458
Fs(q) = m  *  ϕ*A/(m*c)*(1 + q)
display(Fs(0))
# doubles
display(Fs(1))

a_e = 1.49597870700e11
irr0 = 1366
# irr(D) = a_e^2*irr0/D^2;
L(D) = 4*pi*irr0*D^2
display(L(a_e))

# include("../Tools/Newton.jl")
# L2(D) = 4*pi*irr0*D^2 - 3.9e26
# df = newton(a_e, L2, max_iter=1000, eps_N=1e11, eps_jac=1e5)
# display(df)
# display(L(df))

d = sqrt(3.9e26/(4*pi*irr0))
display(d)
display(L(d))

# deceleration at perigee making it circular

