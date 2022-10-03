addpath(genpath(fileparts(which(mfilename))));

K = types;

EA = EulerRot;
PRV = PRVRot

dcm = K.DCM(eye(3,3), 'B0', 'N0');
ea = K.EulerAng([1,pi/4,1],[3,2,1],'B0', 'N0');
prv = K.PRV(pi/2,[1,0,0],'B0', 'N0')

EA.euler2dcm(ea)
EA.dcm2euler(dcm,'ZYX')
EA.euler2om(ea)
EA.om2euler(ea)

PRV.PRV2dcm(prv)
PRV.dcm2PRV(dcm)
PRV.PRV2om(prv)
PRV.om2PRV(prv)