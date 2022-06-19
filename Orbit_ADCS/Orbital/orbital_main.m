
function [tspan, theta, Panel, beta, visibility, Ix, mu] = orbital_main(pointing, n, npoints, T, H, ecc, om, th, M, w, case_adcs, plot_panel, plot_orbit, n_fail)
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
omega_e = (15.04*pi)/(180*3600);% earth rotational speed
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
[r,v] = plotorbit([r0 v0],0,Tend,npoints,mu_e,1,0);

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

theta_g = (greenwich0);     % rotation om_e.*tspan already implemented in the IGRF function
XYZ = zeros(npoints,3);
B_N = zeros(npoints,3);
for k = 1:npoints
    [g, h]   = IGRF13();
    B_N(k,:) = earthmagfield13(r(k,:)', tspan(k,:), g, h, theta_g, 13);
    XYZ(k,:) = B_N(k,:)/norm(B_N(k,:));
end

theta = linspace(0, n*2*pi, npoints);

[Panel, BBx, BBy, BBz] =  pointing(a_e, irr0, npoints, tspan, R, visibility, r, v, XYZ, th, w);

Panel = fail_mode(Panel, n_fail);

if plot_panel == true
    figure
    plot(tspan, BBx)
    xlim([0,max(tspan)])

    figure
    plot(tspan, BBy)
    xlim([0,max(tspan)])

    figure
    plot(tspan, BBz)
    xlim([0,max(tspan)])
T = T*60;

    figure
    plot(tspan/T,visibility,'k','LineWidth',1.5)
    yticks([0 1])
    ylim([-0.25 1.25])
    yticklabels({'no','yes'})
    xl = xlabel('$n_{orbits}$', 'interpreter', 'latex','fontsize',18);
    yl = ylabel('Light', 'interpreter', 'latex','fontsize',18);
    set(get(gca,'ylabel'),'Rotation',0)
    set(yl, 'Units', 'Normalized', 'Position', [-0.1, 0.91, 0]);
    set(xl, 'Units', 'Normalized', 'Position', [0.85, -0.075, 0]);
    box on
    set(gca,'FontSize',18)

    figure
    plot(tspan/T, Panel)
    xlim([0,max(tspan)/T])
    xl = xlabel('$n_{orbits}$', 'interpreter', 'latex','fontsize',18);
    yl = ylabel('G $\left[\mathrm{\frac{W}{m^2}}\right]$', 'interpreter', 'latex','fontsize',18);
    set(get(gca,'ylabel'),'Rotation',0)
    set(yl, 'Units', 'Normalized', 'Position', [-0.1, 0.91, 0]);
    set(xl, 'Units', 'Normalized', 'Position', [0.85, -0.075, 0]);
    legend('Panel 1','Panel 2','Panel 3','Panel 4')
    box on
    set(gca,'FontSize',18)
end

if plot_orbit == true
    figure
    orbit_earth = plot3(r_e(:,1),r_e(:,2),r_e(:,3),'Color','r','LineWidth',2);
    hold on
    orbit_spacecraft = plot3(R(:,1),R(:,2),R(:,3),'Color','b','LineWidth',2);

    % sun = plot3(0,0,0,'o','MarkerFaceColor','k','MarkerEdgeColor','k');
    earth = plot3(nan,nan,nan,'o','MarkerFaceColor','k','MarkerEdgeColor','k');
    sc = plot3(nan,nan,nan,'o','MarkerFaceColor','k','MarkerEdgeColor','k');
    link = plot3(nan,nan,nan,'r','LineWidth',2);

    axis equal

    %%
    for k = 1:round(npoints/1000):npoints
        set(earth,'XDATA',r_e(k,1))
        set(earth,'YDATA',r_e(k,2))
        set(earth,'ZDATA',r_e(k,3))
        
        set(sc,'XDATA',R(k,1))
        set(sc,'YDATA',R(k,2))
        set(sc,'ZDATA',R(k,3))
        
    %     if visibility(k)
    %         set(link,'XDATA',[0 R(k,1)])
    %         set(link,'YDATA',[0 R(k,2)])
    %         set(link,'ZDATA',[0 R(k,3)])
    %     else
    %         set(link,'XDATA',nan)
    %         set(link,'YDATA',nan)
    %         set(link,'ZDATA',nan)
    %     end
        
        drawnow
    end
end

N  = 244;
A1 = 1005*1e-6;
A2 = 11^2*1e-6;
J  = 3.131*1e-4;
dt = tspan(2) - tspan(1);
[Ix, mu]  = mag_2(dt, N, A1, J, BBx, B_N');
toc
end