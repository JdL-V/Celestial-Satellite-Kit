close all
clear
clc

O = [0 0 0]';

i = [1 0 0]';
j = [0 1 0]';
k = [0 0 1]';

A = 1.578e-4;
B = 3.131e-4;
C = 3.093e-4;

I0 = [A 0 0;
      0 B 0;
      0 0 C];

%% ode
npoints = 10000;
tspan = linspace(0, 95*60*1, npoints);
y0 = [i; j; k; 0.01 + rand*eps; rand*eps; rand*eps];

tic
options = odeset('RelTol',1e-13, 'AbsTol',1e-14);
[time, Y] = ode113(@(t,y) FF(t, y, I0), tspan, y0, options);
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
[az, el] = view;
view([az+180 el])
text(i,j,k, {'i','j','k'})
hold on

filename    = 'POCKET_CUBE.STL';
[p,t,tnorm] = import_stl_fast2(filename,1);
p = p/1e2; % scaling
p = [p(:,1) - max(p(:,1))/2,  p(:,2) - max(p(:,2))/2,  p(:,3) - max(p(:,3))/2]; % centering
p2 = p;

% plot_speed = round(npoints/tspan(end)/norm(y0(end-2:end))*0.06);
plot_speed = round(npoints/2000);
for kk = [1:plot_speed:length(tspan) npoints]

I = Y(kk, 1:3)';
J = Y(kk, 4:6)';
K = Y(kk, 7:9)';

for kn = 1:size(p,1)
    p2(kn,:) = ([I J K]'*p(kn,:)')';
end

try
    delete(ppp)
    delete(hQuiver)
catch
end

ppp = patch('vertices',p2,'faces',t,'FaceColor','c','edgecolor','k');
hQuiver = quiver3(O,O,O,I,J,K,'b');
drawnow
title(['t = ' num2str(round(time(kk))) ' s'])
end

function du = FF(~, u, I0)
i = [1 0 0]';
j = [0 1 0]';
k = [0 0 1]';

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
%     t
end