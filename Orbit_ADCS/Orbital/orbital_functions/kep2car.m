function [rr,vv] = kep2car(a,e,i,OM,om,th,mu)

% kep2car.m - Conversion from Keplerian elements to Cartesian coordinates
%
% PROTOTYPE:
% [r,v]=kep2car(a,e,i,OM,om,th,mu)
%
% DESCRIPTION:
% Conversion from Keplerian elements t oCartesian coordinates .Angles in 
% radians.
%
% INPUT:
% a        [1x1]       Semi-major axis           [km]
% e        [1x1]       Eccentricity              [-]
% i        [1x1]       Inclination               [rad]
% OM       [1x1]       RAAN                      [rad]
% om       [1x1]       Pericentre anomaly        [rad]
% th       [1x1]       True anomaly              [rad]
% mu       [1x1]       Gravitational parameter   [km^3/s^2]
%
% OUTPUT:
% r        [3x1]       Position vector           [km]
% v        [3x1]       Velocity vector           [km/s]

% Step 1:

p=a*(1-e^2);
r=p/(1+e*cos(th));

% Step 2:

r_e=r*cos(th);
r_p=r*sin(th);
r_perif=[r_e;r_p;0];

v_e=-sqrt(mu/p)*sin(th);
v_p=sqrt(mu/p)*(e+cos(th));
v_perif=[v_e;v_p;0];

% Step 3:

R1=[1 0 0;0 cos(i) sin(i);0 -sin(i) cos(i)];
R2=[cos(om) sin(om) 0;-sin(om) cos(om) 0;0 0 1];
R3=[cos(OM) sin(OM) 0;-sin(OM) cos(OM) 0;0 0 1];

% Step 4:

R=R2*R1*R3;

% Step 5:

rr=R\r_perif;
vv=R\v_perif;

end

