format compact
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
I2 = [98. 0. 0.; 
      0. 102. 0.; 
      0. 0. 150.];
      
w = [0.1, 0.02, 0.5]';
H = I2*w;
beta1 = acosd(dot(H./norm(H),w./norm(w)))
alfa1 = acosd(dot(w./norm(w),[0,0,1]))

optiondop = rdpset('RelTol',1e-7,'AbsTol',1e-7,'Refine',10);
warning off
u0 = [1., 0., 0., 0., w(1), w(2), w(3)]';
tic
    [tout, uout] = EPdiff(@fF4, linspace(0,60,1e3), u0, optiondop);
toc
warning on

figure
plot(tout, uout(:,1:4))
figure
plot(tout, uout(:,5:7))

beta2 = zeros(1000,1);
alfa2 = zeros(1000,1);
for i = 1:1000
    w = uout(i,5:7)';
    H = I2*w;
    beta2(i) = acosd(dot(H./norm(H),w./norm(w)));
    alfa2(i) = acosd(dot(w./norm(w),[0,0,1]));
end
figure
plot(tout,[beta2, alfa2])


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