format compact
%1
%a
% kinetic energy from photons
%b
% phi*A/mc*(1 + q)
phi = 1e3;
A = 1;
m = 10;
c = 299792458;
Fs = @(q) m  *  phi*A/(m*c)*(1 + q);
Fs(0)
%c
% doubles
Fs(1)

%d
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

%2
% deceleration at perigee making it circular lowering the apogee

%3
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

L2 = eznewton(-1e6, L123x)

L1 = eznewton(1e6, L123x)

L3 = eznewton(r1*1.01, L123x)

%4
%a
T = 86164;
mu = astroConstants(13);
w = 2*pi/T;
Rgeo = (mu/w^2)^(1/3)
% Harmonics








% Re = astroConstants(23)
% a = 7721
% T = 2*pi*sqrt(a^3/mu)

% vmin = sqrt(mu*(2/(Re+600) - 1/a))
% vmax = sqrt(mu*(2/(Re+1600) - 1/a))


% % dcm = DCM([0.9551 0.0052 0.2962; 0.1504 0.8529 -0.5; -0.2552 0.5221 0.8137],'a','b')
% % dcm2euler(dcm,'XYZ')

% eta = @(kp) 1/sqrt(kp)
% f = @(kp) exp(-pi*eta(kp)/sqrt(1-eta(kp)^2))*100 - 4.31

% eznewton(2,f,4000)