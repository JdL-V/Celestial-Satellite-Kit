function bhbs(h, rho_function, a, b, l, v, J, x_cm, M, cd, x_cp_g, M_mag, dip_res, PS, c, q, i)
%% DATOS:
alt = 500e3; %km
R_T = 6378e3; %km
R = alt+R_T; %km
mu = 398.6e12; %m^3/s^2
%rho_0 = 3.725e-12; %kg/m^3
%rho_amb = rho_0*exp(-alt/58515); %kg/m^3
%[Temperature, aaa, Pressure, rho] = atmosisa();
%rho_amb = 7e-11;%kg/m^3 (300km)
rho_amb = 6.967e-13;%kg/m^3 (500km)

v = sqrt(mu/R); %m/s (Velocidad en órbita circular)
T = 2*pi*R^(3/2)/sqrt(mu); %s (Periodo en órbita circular)

a = 50e-3; %m
b = 50e-3; %m
L = 81e-3; %m

I_x = 380.1e-6; %kg*m^2
I_y = 367e-6; %kg*m^2
I_z = 203.8e-6; %kg*m^2

x_cm = 0.01; %m (Desplazamiento del centro de masas)

M_pl = 0.080; %kg
M_tot = 0.374; %kg

% Desviación máxima del eje z desde la vertical local (coordenadas cuerpo y coordenadas navegación local):
theta = 1*pi/180; %rad (modo normal)
%theta = 30*pi/180; %rad (modo opcional objetivo de oportunidad)

% BOOM (acero):
M_b = M_pl - 0.040; %kg
Radius = ((M_b/7850)/(4/3*pi))^(1/3) %m
I_b = 2/5*M_b*Radius^2;
C_d_b = 2.2; %(Coeficiente de arrastre a 300km)
A_b = 0*pi*Radius^2; %m^2

% SISTEMA (centro de masas del pocketqube en centro geometrico)
%d = 0.200; %longitud del boom
%d = 0.300; %longitud del boom
d = 0.250; %longitud del boom
x_cg_s = (0*(M_tot-M_b)+(d+L/2)*M_b)/(M_tot);
%I_x_s = 195.5e-6;
%I_y_s = 2e-3;
%I_z_s = 2e-3;

 %I_x_s = 195.5e-6;
 %I_y_s = 4e-3;
 %I_z_s = 4e-3;

I_x_s = 195.5e-6;
I_y_s = 3e-3;
I_z_s = 3e-3;


%% MOMENTO POR GRADIENTE GRAVITATORIO:
Mg = 3*mu/(2*R^3)*abs(I_z-I_x)*sin(2*theta) %N*m


%% MOMENTO AERODINÁMICO:
C_d = 2.2; % (Coeficiente de arrastre a 300km)
x_cp_g = 0.07; % (Centro de presiones con respecto a la gravedad)
A = a*L; %m^2 (Área frontal max)
%A = a*b; %m^2 (Área frontal min)

Fa = 0.5*(rho_amb*C_d*A*v^2); %N

Ma = Fa*x_cp_g %N*m


%% MOMENTO MAGNÉTICO:
Mm_T = 7.96e15; %T*m^3
dip_res = 0.01; %A*m^2 (Dipolo residual a 400km)
B_T = 2*Mm_T/R^3; %T (Aproximación para una órbita polar )

Mm = dip_res*B_T %N*m



%% MOMENTO POR RADIACIÓN SOLAR:
Ps = 1.371; %W/m^2 (Constante solar)
c = 3e8; %m/s (Velocidad de la luz)
q = 0.6; % (Coeficiente de reflexión [0,1])
i = 0*pi/180; %rad (Ángulo de incidencia del Sol)
%x_cps_g = ; %m (Centro de presión solar con respecto al de gravedad)

%Fs = Ps/c*(1+q)*A*cos(i); %N

%Ms = Fs*(x_cps_g) %N*m


%% MOMENTOS DEL SISTEMA:
Mg_s = 3*mu/(2*R^3)*abs(I_z_s-I_x_s)*sin(2*theta) %N*m

%Fa = 0.5*(rho_amb*C_d*A*v^2); %N
Fa_b = 0.5*(rho_amb*C_d_b*A_b*v^2); %N

Ma_s = abs(Fa*x_cg_s-Fa_b*(d+L/2-x_cg_s)) %N*m
