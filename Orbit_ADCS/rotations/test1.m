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
b = [2  2  2  2, -2 -2 -2 -2, 2  2 -2 -2,  2  2 -2 -2,  2  2 -2 -2,  2  2 -2 -2;
     1  1 -1 -1,  1  1 -1 -1, 1  1  1  1, -1 -1 -1 -1,  1 -1 -1  1,  1 -1 -1  1;
     1 -1 -1  1,  1 -1 -1  1, 1 -1 -1  1,  1 -1 -1  1,  1  1  1  1, -1 -1 -1 -1]*0.2;
a = zeros(size(b));

for kk = [1:round(npoints/1000):length(tspan) npoints]
    %1:round(npoints/tspan(end)/norm(y0(end-2:end))*0.06):length(tspan)

I = Y(kk, 1:3)';
J = Y(kk, 4:6)';
K = Y(kk, 7:9)';

for kn = 1:24
    a(:,kn) = [I J K]'*b(:,kn);
end
try
    delete(ppp1); delete(ppp2); delete(ppp3); delete(ppp4); delete(ppp5); delete(ppp6)
    delete(hQuiver)
catch
end
ppp1 = patch(a(1, 1:4), a(2, 1:4), a(3, 1:4), 'c');
ppp2 = patch(a(1, 5:8), a(2, 5:8), a(3, 5:8), 'c');
ppp3 = patch(a(1, 9:12), a(2, 9:12), a(3, 9:12), 'c');
ppp4 = patch(a(1, 13:16), a(2, 13:16), a(3, 13:16), 'c');
ppp5 = patch(a(1, 17:20), a(2, 17:20), a(3, 17:20), 'c');
ppp6 = patch(a(1, 21:24), a(2, 21:24), a(3, 21:24), 'c');

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
    
    M = [1e-4*I'*i J'*j+K'*j J'*k+K'*k]'*1e-5;

    du = [dI; dJ; dK; dw + M];
%     t
end