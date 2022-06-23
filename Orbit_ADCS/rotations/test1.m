close all
clear
clc

O = [0 0 0]';

i = [1 0 0]';
j = [0 1 0]';
k = [0 0 1]';

A = 1;
B = 2;
C = 2;

I0 = [A 0 0;
      0 B 0;
      0 0 C];

tspan = linspace(0, 3600, 10000);
y0 = [i; j; k; 1; 1; 0];

tic
options = odeset('RelTol',1e-13, 'AbsTol',1e-14);
[~, Y] = ode113(@(t,y) FF(t, y, I0), tspan, y0, options);
toc
%% plot
figure
quiver3(O,O,O,i,j,k,'k')
xlabel('x')
ylabel('y')
zlabel('z')
axis equal
xlim([-1 1]);
ylim([-1 1]);
zlim([-1 1]);
hold on
for kk = 1:length(tspan)
% om = deg2rad(1);

% Rz = [1 0 0
% 0 cos(om*t(kk)) -sin(om*t(kk));
% 0 sin(om*t(kk)) cos(om*t(kk))];

% I = Rz*i;
% J = Rz*j;
% K = Rz*k;
I = Y(kk, 1:3)';
J = Y(kk, 4:6)';
K = Y(kk, 7:9)';
hQuiver = quiver3(O,O,O,I,J,K);
drawnow
pause(0.05)
delete(hQuiver)
end

function du = FF(t, u, I0)
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
    
    du = [dI; dJ; dK; dw];
    t
end