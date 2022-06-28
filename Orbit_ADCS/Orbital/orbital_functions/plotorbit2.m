function [r,v] = plotorbit2(I0, solver, y0, tstart, tfinal, npoints, mu, flag, plot, varargin)

% plotorbit.m - propagate the orbit from the initial state y0 in the timespan from tstart to tfinal.
%
% PROTOTYPE:
% [r,v,orbit]=plotorbit(y0,tstart,tfinal,npoints,mu,flag,varargin).
%
% DESCRIPTION:
% Function to propagate either non pertubed or pertubed(J2 effect)
% orbits.Two body problem. NOTE: J2 effect implemented only if the central body is the earth
%
% INPUT:
% y0        [1x6]   vector that contains the initial position and the initial
%                   velocity.Where r0 [1x3] and v0 [1x3] then y0=[r0,v0].
% tstart    [1x1]   starting time point of integration.
% tfinal    [1x1]   final time point of integration.
% npoints   [1x1]   number of point used to integrate the ode.
% mu        [1x1]   gravitational parameter of the central body
% flag      [1x1]   must be 1 or 0. 1 to propagate the orbit accounting J2
%                   effect (valid only if the central mass attractor is the
%                   Earth, 0 to propagate the orbit without any disturbance.
% plot      [1x1]   If 1 it plots the orbit,if 0 it doesn't plot it.
% varargin  [1x?]   contains the desired property and corresponding values
%                   for the line.
%
% OUTPUT:
% r     [npointsx3] vector that contains all the position in the timespan
% v     [npointsx3] vector that contains all the velocity in the timespan
% orbit [1x1 line]  line containing the orbit.(used for further operations)
%
% EXAMPLE:
% mu=398600;
% y0=[-1.829361896546810e+04,3.189584445999387;-1.703199789951243e+03,-3.112192033494531;-6.658332780438067e+03,1.160913798022232];
% tstart=0;
% tfinal=5.778137031949276e+4;
% npoints=500;
% [r,v,orbit]=plotorbit(y0,tstart,tfinal,npoints,mu,1,'Color','red','LineWidth',2.5,'Marker','o','MarkerSize',7,'MarkerFaceColor','g');
% grid on
% axis equal
%
%   Notes for upgrading this function:
%       It is possible to change some operations to make it faster.
%       - DO NOT change the structure of the function, as well as its
%           prototype.
%       - DO NOT change the output type.
%       - DO NOT change the input to be given in degrees
%       Contact the authors for modification: jmjorgem11@gmail.com
%
% CALLED FUNCTIONS:
% keplerian_orbit,pertubed_keplerian_orbit(added at the bottom of this file)
%
% CHANGELOG:
%  23/11/2020, Riccardo Majer, Jorge Martínez: creation.
%  01/06/2022, Jorge Martínez: update - commented displays.
%  01/06/2022, Jorge Martínez: update - norm computed once (r_s) and fixed J2.
% ------------------------------------------------------------------------------------------------------------

switch flag
    case 1
        % for orbit propagation accounting J2 perturbation. NOTE:J_2 and R_e in the if section are valid only around the earth
    case 0
        % for orbit propagation without disturbances
    otherwise
        error('flag must be 1 or 0')
end

options=odeset('RelTol',1e-13,'AbsTol',1e-14);

tspan=linspace(tstart,tfinal,npoints);
if flag
    R_e=6378.137;
    J_2=0.00108263;
    [~,Y]=solver(@(t,y) pertubed_keplerian_orbit(t,y,mu,J_2,R_e,I0),tspan,y0,options);
%     disp('pertubed');
else
    [~,Y]=solver(@(t,y) keplerian_orbit(t,y,mu),tspan,y0,options);
%     disp('normal');
end
r=[Y(:,1),Y(:,2),Y(:,3)];
v=[Y(:,4),Y(:,5),Y(:,6)];
if plot
    plot3(r(:,1),r(:,2),r(:,3),varargin{:});
end
end

%CALLED FUNCTIONS

function dy=keplerian_orbit(~,y,mu)
% keplerian_orbit ODE system for the resolution of the two body problem
%
% PROTOTYPE:
% dy=keplerian_orbit(t,y,mu)
%
% INPUT:
% t  [1x1] [T] Time(can be omitted, as the system is autonomous)
% y  [1x6] [L,L/T] state of s/c(position and velocity[three axis])
% mu [1x1] gravitational parameter of primary body
%
% OUTPUT:
% dy [1x6] derivative of the state
r_s = norm([y(1), y(2), y(3)]);

dy=[y(4) % r_x
    y(5) % r_y
    y(6) % r_z
    -mu/r_s^3.*y(1) % v_x
    -mu/r_s^3.*y(2) % v_y
    -mu/r_s^3.*y(3) % v_z
    ];
end

function dy=pertubed_keplerian_orbit(~,y,mu,J_2,R_e,I0)
% pertubed_keplerian_orbit ODE system for the resolution of the pertubed two body problem
%
% PROTOTYPE:
% dy=perubed_keplerian_orbit(t,y,mu,J_2,R_e)
%
% INPUT:
% t   [1x1] [T] Time(can be omitted, as the system is autonomous)
% y   [1x6] [L,L/T] state of s/c(position and velocity[three axis])
% mu  [1x1] gravitational parameter of primary body
% J_2 [1x1]
% R_e [1x1] [L] Earth's radius
%
% OUTPUT:
% dy [1x18] derivative of the state
r_s = norm([y(1), y(2), y(3)]);
A_J2 = 0.5*J_2*(R_e/r_s)^2;

dy=[y(4) % r_x
    y(5) % r_y
    y(6) % r_z
    -mu/r_s^3*y(1)*(1 + 3*A_J2*(1 - 5*(y(3)/r_s)^2))   % v_x
    -mu/r_s^3*y(2)*(1 + 3*A_J2*(1 - 5*(y(3)/r_s)^2))   % v_y
    -mu/r_s^3*y(3)*(1 + 3*A_J2*(3 - 5*(y(3)/r_s)^2))]; % v_z

u = y(7:18)

% rotations
I = u(1:3);
J = u(4:6);
K = u(7:9);
w = u(10:12);

A = I0(1,1);
B = I0(2,2);
C = I0(3,3);

dI = cross(w, I);
dJ = cross(w, J);
dK = cross(w, K);
    
dw = -[(C - B)/A*w(2)*w(3);
          (A - C)/B*w(1)*w(3);
          (B - A)/C*w(2)*w(1)];
    
M = [1e-4*I'*i J'*j+K'*j J'*k+K'*k]'*1e-8;

du = [dI; dJ; dK; dw + M./[A B C]'];

dy = [dy du'];
end