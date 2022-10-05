addpath(genpath(fileparts(which(mfilename))));
clc
clear
close all

TP = types;
global ROT
ROT = RotModule;

dcm = TP.DCM(eye(3,3),'B0', 'N0');
ea  = TP.EulerAng([1,pi/4,1],[3,2,1],'B0', 'N0');
prv = TP.PRV(pi/2,[1,0,0],'B0', 'N0');
ep  = TP.quaternion([0.6301, 1.0502, -0.6301, -1.0502],'B0', 'N0');
crp = TP.CRP([0.1, 0.2, 0.3],'B0', 'N0');
mrp = TP.MRP([0.1, 0.2, 0.3],'B0', 'N0');

ROT.euler2dcm(ea);
ROT.dcm2euler(dcm,'ZYX');
ROT.euler2om(ea);
ROT.om2euler(ea);

ROT.PRV2dcm(prv);
ROT.dcm2PRV(dcm);
ROT.PRV2om(prv);
ROT.om2PRV(prv);

ROT.EP2dcm(ep);
ROT.dcm2EP(dcm, 0);
ROT.EP2PRV(ep);
ROT.PRV2EP(prv);
ROT.SumEP(ep,ep);
ROT.EP2om(ep.x);
ROT.om2EP(ep.x);

ROT.CRP2dcm(crp);
ROT.dcm2CRP(dcm);
ROT.CRP2PRV(crp);
ROT.PRV2CRP(prv);
ROT.SumCRP(crp,crp);
ROT.CRP2om(crp.x);
ROT.om2CRP(crp.x) ;

ROT.MRP2dcm(mrp);
ROT.dcm2MRP(dcm);
ROT.MRP2PRV(mrp);
ROT.PRV2MRP(prv);
ROT.MRP2CRP(mrp);
ROT.CRP2MRP(crp);
ROT.SumMRP(mrp,mrp);
ROT.MRP2om(mrp.x);
ROT.om2MRP(mrp.x);

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

optiondop = rdpset('RelTol',1e-7,'AbsTol',1e-7,'Refine',10);

[tout, uout] = ROT.EPdiff(@fq2,linspace(0,60,1e3),[1 0 0 0],optiondop);
% a = [];
% for nnn = 1:size(uout,1)
%     a = [a norm(uout(nnn,:)) - 1];
% end
% max(abs(a))
plot(tout,uout)

Om0 = 0.5;
u0 = [0., 0., 0., 3/5*Om0, 0., 4/5*Om0];

[tout, uout] = ROT.MRPdiff(@fF3,linspace(0,60,1e3),u0,optiondop);
figure
plot(tout,uout(:,1:3))
figure
plot(tout,uout(:,4:6))

function var = fq2(t, u)
    global ROT
    w = @(t) deg2rad(50).*[sin(0.1*t), 0.01, cos(0.1*t)]';
    M = ROT.om2EP(u);
    var = M*w(t);
end

function var = fF3(t, u)
    global ROT
    m = 1.;
    R = 1.;
    Ib = [4*m*R^2/3 + m*R^2/4 0 0; 0 4*m*R^2/3 + m*R^2/4 0; 0 0 m*R^2/2];
    
    Lc = [0,0,0]';
    var = [ROT.om2MRP(u(1:3))*u(4:6); -cross2mat(u(4:6))*Ib*u(4:6) + Lc];
end