format compact
close all
%1
m = 1.;
R = 1.;
Om0 = 0.5;
% H = 4R
Ib = [4*m*R^2/3 + m*R^2/4   0                       0;
      0                     4*m*R^2/3 + m*R^2/4     0;
      0                     0                       m*R^2/2];

%2

global I2
I2 = [100. 0. 0.; 
      0. 100. 0.; 
      0. 0. 150.];
      
w = [0.1, 0.02, 0.5]';
H = I2*w;
ht = norm(H(1:2))
h = norm(H)
nut_ang = asind(ht/h)
prec_rate = rad2deg(h/I2(1,1))

%3
%Polhode plot

%4
I2 = [210. 0. 0.; 
      0. 200. 0.; 
      0. 0. 118.];
w = [0.05, 0.02, -0.02]';

warning off
u0 = [0., 0., 0., w(1), w(2), w(3)]';
tic
    [tout, uout] = dop853(@fF5, linspace(0,100,1e3), u0, optiondop);
toc

figure
plot(tout, uout(:,1:3))
figure
plot(tout, uout(:,4:6))






% I2 = [20e6 0.   0.; 
%       0.   20e6 0.; 
%       0.   0.   5e6];
      
% w3 = [0.0262]';
% h3 = I2(3,3)*w3;

% f = @(x) sind(30) - x/(sqrt(x^2 + h3^2))
% ht = eznewton(w3, f, 10000)
% h = sqrt(ht^2 + h3^2)
% prec_rate = rad2deg(h/I2(1,1))
% spin = rad2deg(prec_rate*(20e6 - 5e6)*cosd(30)/5e6)

% wt = prec_rate*sind(30)
% ww = norm([wt wt w3])
% w = [wt 0 w3]
% H = [ht 0 h3]
% e3 = [0 0 1]
% acosd(dot(w,e3)/norm(w)/norm(1))



function du = fF4(t,u)
    global I2
    Lc = [0,0,0]';
    du = [om2EP(u(1:4))*u(5:7); I2\(-skewsym(u(5:7))*I2*u(5:7) + Lc)];
end

function du = fF5(t,u)
    global I2
    Lc = [0,0,0]';
    du = [om2euler(EulerAng(u(1:3), [3,2,1], "A", "A"))*u(4:6); I2\(-skewsym(u(4:6))*I2*u(4:6) + Lc)];
end

