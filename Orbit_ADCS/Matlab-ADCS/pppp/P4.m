M = 1e3;
At = 1e4;
ps = 1361/299792458;
as = 2*ps*At/M;
f = @(t) as/2*t^2 - 4e8;
newton(1e6, f, 1000, 1e-6, 1e-8)
eznewton(1e6, f, 1000)
tt = sqrt(4e8*2/as)

%%%%%%%%%%%%%%%%%%%%%

I3 = 100;
I1 = 250;
I2 = I1;
global b
global Ic
b = 5;

Ic = [I1 0 0; 0 I2 0; 0 0 I3];
u0 = [8 9 7];
warning off
optiondop = rdpset('RelTol',1e-14,'AbsTol',1e-14,'Refine',10);
[t, u] = dop853(@fF3,linspace(0,60,1e5),u0,optiondop);
plot(t,u)

% dw3 = -b*w3/I33
% w3 = K*e^(-b*t/I33)

function du = fF3(t, u)
    global b
    global Ic
    Lc = -b.*u;
    du = Ic\(-skewsym(u)*Ic*u + Lc);
end