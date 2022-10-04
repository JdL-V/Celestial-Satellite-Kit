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