
function [tspan, theta, Panel, beta, visibility, Ix, mu] = orbital_main(pointing, n, npoints, T, H, ecc, om, th, M, w, plot_panel, plot_sunorbit, plot_earthorbit, n_fail)
tic
%% Costants and initial conditions:
AU = astroConstants(2);         % distance 1AU in km

th_e = deg2rad((M + 1.5)*360/12);       % theta earth
irr0 = 1366;
% irr = irr0 + 45*cos((M + 1.5)*2*pi/12) % earth irradiance model

% sun centered body:
mu_s = astroConstants(4);
position_s = [0;0;0];

% earth centered body:
mu_e = astroConstants(13);      % mu earth
Re_e = astroConstants(23);      % R earth
% omega_e = (15.04*pi)/(180*3600);% earth rotational speed
greenwich0 = 0;                 % initial Greenwich position

% earth orbit around sun:
a_e = 1*AU;
ecc_e = 0.01671022;
inc_e = deg2rad(0.00005);
raan_e = deg2rad(-11.26064);
om_e = deg2rad(102.94719);

% simulation time and number of points:
Tend = n*T*60;
tspan = linspace(0,Tend,npoints)';

%% CALCULATIONS:

% earth initial position around sun:
[r0_e,v0_e] = kep2car(a_e,ecc_e,inc_e,raan_e,om_e,th_e,mu_s);

% S/C initial position around earth:
a = (mu_e*(T*60/(2*pi))^2)^(1/3);

% SSO case
J2 = astroConstants(9);
p = a;
d_om_sun = 2*pi/(365.2411984*24*3600);
inc = acos(d_om_sun/(-3/2*J2*Re_e^2/p^2*sqrt(mu_e/a^3)));    % orbital plane inclination SSO.
raan_0 = deg2rad(((H - 12)*90/6) - 11.26064);                  % raan_e corrected
raan = raan_0 + d_om_sun*(M - 1)*(365.2411984*24*3600)/12;   % orbital RAAN SSO

[r0,v0] = kep2car(a,ecc,inc,raan,om,th,mu_e);

% earth orbit propagation around sun:
[r_e,~] = plotorbit([r0_e v0_e],0,Tend,npoints,mu_s,0,0);

% S/C orbit propagation around earth:
[r,v] = plotorbit([r0 v0],0,Tend,npoints,mu_e,0,0);

% S/C orbit propagation around sun:
R = r + r_e;

%% Eclipses
d_earth = zeros(npoints,1); % minimum distance vector sun-S/C segment

for k = 1:npoints
    d_earth(k) = point_to_line([r_e(k,1);r_e(k,2);r_e(k,3)], position_s, [R(k,1);R(k,2);R(k,3)]);
end

visibility = logical(d_earth > Re_e);

%% Beta angle
beta = zeros(npoints,1); % 
for k = 1:npoints
    sun_dir = (R(k,:)/norm(R(k,:)))';
    ht = cross(r(k,:), v(k,:));
    hn = ht*sun_dir/(norm(ht)*norm(sun_dir));
    beta = rad2deg(acos(hn));
    if beta > 90
        beta = beta - 90;
    end
end

%% Magnetic field:
theta_g = (greenwich0);     % rotation omega_e.*tspan already implemented in the IGRF function
XYZ = zeros(npoints,3);
B_N = zeros(npoints,3);
for k = 1:npoints
    [g, h]   = IGRF13();
    B_N(k,:) = earthmagfield13(r(k,:)', tspan(k,:), g, h, theta_g, 13);
    XYZ(k,:) = B_N(k,:)/norm(B_N(k,:));
end

theta = linspace(0, n*2*pi, npoints);

[Panel, BBx, BBy, BBz] = pointing(a_e, irr0, npoints, tspan, R, visibility, r, v, XYZ, th, w);

Panel = fail_mode(Panel, n_fail);

if plot_panel == true
    panel_plot(BBx, BBy, BBz, tspan, visibility, T, Panel)
end

if plot_earthorbit == true
    earth_sat(a, r, tspan, npoints, 0.6)
end

if plot_sunorbit == true
    sun_orbit(r_e, R, npoints)
end

N  = 244;
A1 = 1005*1e-6;
A2 = 11^2*1e-6;
J  = 3.131*1e-4;
dt = tspan(2) - tspan(1);
[Ix, mu]  = mag_2(dt, N, A1, J, BBx, B_N');

figure
plot(tspan, [mu(1) mu mu(end)])
title('mag2')

N = [250; 850; 850]';
A = [33^2; 11^2; 11^2]'*1e-6;
I = [1; 1; 1]'*2e-3;

[T_mx, T_my, T_mz] = mag_3(dt, tspan, J, BBx, BBy, BBz, B_N', w, N, A, I);
toc

% blue: x // orange: y // yellow: z
figure
plot(tspan, T_mx);
title('mag3 x')

figure
plot(tspan, T_my);
title('mag3 y')

figure
plot(tspan, T_mz);
title('mag3 z')

% figure
% plot(tspan, sqrt(sum(T_mx.^2, 2)/3));
% title('mag3 x sqrt')

% figure
% plot(tspan, sqrt(sum(T_my.^2, 2)/3));
% title('mag3 y sqrt')

% figure
% plot(tspan, sqrt(sum(T_mz.^2, 2)/3));
% title('mag3 z sqrt')

end