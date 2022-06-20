clear
close all
clc
format compact

addpath(genpath(fileparts(which(mfilename))));


%***************************************************************************************************************************************
%***************************************************************************************************************************************
%**               VARIABLES A MODIFICAR:                                                                                              **
%***************************************************************************************************************************************
%***************************************************************************************************************************************

n = 200;                   % number of orbits (15)
dT = 50;                   % time step (0.5)
case_adcs = 1;            % case ADCS (1=sigue las lineas de campo, 2=x mirando a nadir, 3=x perpendicular al plano, 4=paneles al sol)        
ciclo = 1;                % ciclo de consumos de los componentes (1 = 30s, 2 = 60s) 
n_fail = [];               % failed panel vector (none = [], 1+3+4 = [1 3 4]) (zp=1 zn=2 yp=3 yn=4)
H = 12;                  % Hora solar [horas]
w = 0.1;                  % velocidad de rotación [rad/s] (0.01,0.1)

%***************************************************************************************************************************************
%***************************************************************************************************************************************



% spacecraft orbit around earth:
h = 500;
mu_e = astroConstants(13);      % mu earth
Re_e = astroConstants(23);      % R earth
a = Re_e + h;
T = round(2*pi*sqrt(a^3/mu_e)/3)/20;                                 % orbit period S/C around earth [min]
ecc = 0;
om = deg2rad(0);
th = deg2rad(0);
M = 1;                                  % Month of the year March = 1; Feb = 12

% simulation time and number of points:

npoints = n*T*60/dT + 1;

pointing = adcs_select(case_adcs);

[tspan, theta, Panel, beta, visibility, I, mu] = orbital_main(pointing, n, npoints, T, H, ecc, om, th, M, w, true, false, false, n_fail);

% tipo_panel = 2;           % panel type (1 = galio, 2 = silicio (inutil)) 
% eta_regulador = 0.5;      % rendimiento regulador (0.3-0.6)

% I = [I(1) I I(end)];
% rho = 1;
% N = 250;
% L   = N*116e-3;
% R   = rho*L;
% V = R*I;

% eta_dcdcx = 0.5;
% fprintf('Potencia media = %d [mW]\n', mean(I.*V)*1e3/eta_dcdcx)
% fprintf('Potencia min = %d [mW]\n', min(I.*V)*1e3/eta_dcdcx)
% fprintf('Potencia max = %d [mW]\n', max(I.*V)*1e3/eta_dcdcx)

% figure
% plot(tspan, I)

% figure
% plot(tspan, V)

% figure
% plot(tspan, R*I.^2*1e3/eta_dcdcx)

% mu = [mu(1) mu mu(end)];
% figure
% plot(tspan/T/60, mu,'k','linewidth',1)

% mun = 250*1/R*1005*1e-6;

% xx = trapz(tspan, mu)
% xlim([0,max(tspan)/T/60])
% xl = xlabel('$n_{orbits}$', 'interpreter', 'latex','fontsize',18);
% yl = ylabel('$\tilde{\mu}$ $\left[\mathrm{A \cdot m^2}\right]$', 'interpreter', 'latex','fontsize',18);
% set(get(gca,'ylabel'),'Rotation',0)
% set(yl, 'Units', 'Normalized', 'Position', [-0.075, 0.91, 0]);
% set(xl, 'Units', 'Normalized', 'Position', [0.85, -0.075, 0]);
% box on
% set(gca,'FontSize',18)

% xx/mun/(tspan(end))*100