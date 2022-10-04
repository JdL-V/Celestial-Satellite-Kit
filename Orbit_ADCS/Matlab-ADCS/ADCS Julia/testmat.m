addpath(genpath(fileparts(which(mfilename))));

K = types;

EA = EulerRot;
PRV = PRVRot;
EP = EPRot;
CRP = CRPRot;

dcm = K.DCM(eye(3,3),'B0', 'N0');
ea = K.EulerAng([1,pi/4,1],[3,2,1],'B0', 'N0');
prv = K.PRV(pi/2,[1,0,0],'B0', 'N0');
ep = K.quaternion([0.6301, 1.0502, -0.6301, -1.0502],'B0', 'N0');
crp = K.CRP([0.1, 0.2, 0.3],'B0', 'N0');
mrp = K.MRP([0.1, 0.2, 0.3],'B0', 'N0');

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