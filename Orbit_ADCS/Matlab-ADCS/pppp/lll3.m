format compact
% kinetic energy from photons
% phi*A/mc*(1 + q)
phi = 1e3;
A = 1;
m = 10;
c = 299792458;
Fs = @(q) m  *  phi*A/(m*c)*(1 + q);
Fs(0)
% doubles
Fs(1)

a_e = 1.49597870700e11;
irr0 = 1366;
% irr(D) = a_e^2*irr0/D^2;
L = @(D) 4*pi*irr0*D^2;
L(a_e)

% include("../Tools/Newton.jl")
% L2(D) = 4*pi*irr0*D^2 - 3.9e26
% df = newton(a_e, L2, max_iter=1000, eps_N=1e11, eps_jac=1e5)
% display(df)
% display(L(df))

d = sqrt(3.9e26/(4*pi*irr0))
L(d)

% deceleration at perigee making it circular




r1 = 741.35e6;
r2 = 0.;

G = astroConstants(1);
m1 = 1.989e30;
m2 = 1.898e27;
T = 2*pi*sqrt(r1^3/(G*m1));
om = 2*pi/T;

L123x = @(x) G*m1*(x - r1)/abs(x - r1)^3 + G*m2*(x + r2)/abs(x + r2)^3 - om^2*x;
L4y = (abs(r1) + abs(r2))*sind(60)
L4x = (abs(r1) + abs(r2))*cosd(60)

L2 = newton(-1e6, L123x, 1000, 1e-14, 1e-4)

L1 = newton(1e6, L123x, 1000, 1e-14, 1e-4)

L3 = newton(r1*1.01, L123x, 1000, 1e-14, 1e-4)