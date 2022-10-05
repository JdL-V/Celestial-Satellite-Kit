addpath(genpath(fileparts(which(mfilename))));

TP = types;

EA  = EulerRot;
PRV = PRVRot;
EP  = EPRot;
CRP = CRPRot;
MRP = MRPRot;

dcm = TP.DCM(eye(3,3),'B0', 'N0');
ea  = TP.EulerAng([1,pi/4,1],[3,2,1],'B0', 'N0');
prv = TP.PRV(pi/2,[1,0,0],'B0', 'N0');
ep  = TP.quaternion([0.6301, 1.0502, -0.6301, -1.0502],'B0', 'N0');
crp = TP.CRP([0.1, 0.2, 0.3],'B0', 'N0');
mrp = TP.MRP([0.1, 0.2, 0.3],'B0', 'N0');

EA.euler2dcm(ea)
EA.dcm2euler(dcm,'ZYX')
EA.euler2om(ea)
EA.om2euler(ea)

PRV.PRV2dcm(prv)
PRV.dcm2PRV(dcm)
PRV.PRV2om(prv)
PRV.om2PRV(prv)

EP.EP2dcm(ep)
EP.dcm2EP(dcm, 0)
EP.EP2PRV(ep)
EP.PRV2EP(prv)
EP.SumEP(ep,ep)
EP.EP2om(ep.x)
EP.om2EP(ep.x)

CRP.CRP2dcm(crp)
CRP.dcm2CRP(dcm)
CRP.CRP2PRV(crp)
CRP.PRV2CRP(prv)
CRP.SumCRP(crp,crp)
CRP.CRP2om(crp.x)
CRP.om2CRP(crp.x) 

MRP.MRP2dcm(mrp)
MRP.dcm2MRP(dcm)
MRP.MRP2PRV(mrp)
MRP.PRV2MRP(prv)
MRP.MRP2CRP(mrp)
MRP.CRP2MRP(crp)
MRP.SumMRP(mrp,mrp)
MRP.MRP2om(mrp.x)
MRP.om2MRP(mrp.x)

v1b = [0.8273, 0.5541, -0.0920]'
v2b = [-0.8285, 0.5522, -0.0955]'
v1n = [-0.1517, -0.9669, 0.2050]'
v2n = [-0.8393, 0.4494, -0.3044]'

tic
t1 = TRIAD(v1b, v2b, v1n, v2n)
toc

tic
% t2 = qMethod([v1b v2b], [v1n v2n],[1., 1.])
toc

tic
t3 = QUEST([v1b v2b], [v1n v2n],[1., 1.])
toc

tic
t4 = OLAE([v1b v2b], [v1n v2n],[1., 1.])
toc