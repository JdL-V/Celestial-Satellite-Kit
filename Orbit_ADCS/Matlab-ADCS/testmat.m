clc
clear
close all

t0 = cputime;
addpath(genpath(fileparts(fileparts(which(mfilename)))));

dcm = DCM(eye(3,3),'B0', 'N0');
ea  = EulerAng([1,pi/4,1],[3,2,1],'B0', 'N0');
prv = PRV(pi/2,[1,0,0],'B0', 'N0');
ep  = quat([0.6301, 1.0502, -0.6301, -1.0502],'B0', 'N0');
crp = CRP([0.1, 0.2, 0.3],'B0', 'N0');
mrp = MRP([0.1, 0.2, 0.3],'B0', 'N0');

euler2dcm(ea);
dcm2euler(dcm,'ZYX');
euler2om(ea);
om2euler(ea);

PRV2dcm(prv);
dcm2PRV(dcm);
PRV2om(prv);
om2PRV(prv);

EP2dcm(ep);
dcm2EP(dcm);
sheppard(dcm);
EP2PRV(ep);
PRV2EP(prv);
SumEP(ep,ep);
EP2om(ep.x);
om2EP(ep.x);

CRP2dcm(crp);
dcm2CRP(dcm);
CRP2PRV(crp);
PRV2CRP(prv);
SumCRP(crp,crp);
CRP2om(crp.x);
om2CRP(crp.x);

MRP2dcm(mrp);
dcm2MRP(dcm);
MRP2PRV(mrp);
PRV2MRP(prv);
MRP2CRP(mrp);
CRP2MRP(crp);
SumMRP(mrp,mrp);
MRP2om(mrp.x);
om2MRP(mrp.x);

v1b = [0.8273, 0.5541, -0.0920]';
v2b = [-0.8285, 0.5522, -0.0955]';
v1n = [-0.1517, -0.9669, 0.2050]';
v2n = [-0.8393, 0.4494, -0.3044]';

tic
    t1 = TRIAD(v1b, v2b, v1n, v2n).Mat;
toc

tic
    t2 = qMethod([v1b v2b], [v1n v2n],[1., 1.]).Mat;
toc

tic
    t3 = QUEST([v1b v2b], [v1n v2n],[1., 1.]).Mat;
toc

tic
    t4 = OLAE([v1b v2b], [v1n v2n],[1., 1.]).Mat;
toc
warning off
optiondop = rdpset('RelTol',1e-7,'AbsTol',1e-7,'Refine',10);

tic
    [tout, uout] = EPdiff(@EP_testDF,linspace(0,60,1e2),[1 0 0 0],optiondop);
    % figure; plot(tout,uout)
toc

tic
    Om0 = 0.5;
    u0 = [0., 0., 0., 3/5*Om0, 0., 4/5*Om0];
    [tout, uout] = MRPdiff(@MRP_testDF,linspace(0,60,1e3),u0,optiondop);
    % figure; plot(tout,uout(:,1:3))
    % figure; plot(tout,uout(:,4:6))
toc
warning on

xs = 1:100;
ys = sin(xs);
xout = 1:0.1:100;
warning off
yout = LGinterp(xs, ys, xout, 8);
warning on
% figure
% plot(xs, ys, '.')
% hold on
% plot(xout, yout)
disp(cputime - t0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function var = EP_testDF(t, u)
    global ROT
    w = @(t) deg2rad(50).*[sin(0.1*t), 0.01, cos(0.1*t)]';
    M = om2EP(u);
    var = M*w(t);
end

function var = MRP_testDF(t, u)
    global ROT
    m = 1.;
    R = 1.;
    Ib = [4*m*R^2/3 + m*R^2/4 0 0; 0 4*m*R^2/3 + m*R^2/4 0; 0 0 m*R^2/2];
    
    Lc = [0,0,0]';
    var = [om2MRP(u(1:3))*u(4:6); -skewsym(u(4:6))*Ib*u(4:6) + Lc];
end