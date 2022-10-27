ROT = RotModule;
TP = types;
format compact

% 1
% a Check orthonormal reference frame
a1 = [0.5, sqrt(3)/2, 0]';
a2 = [0., 0., 1.]';
a3 = [sqrt(3)/2, -0.5, 0]';
cross(a1,a2)==a3

b1 = [0, 1, 0]';
b2 = [1, 0, 0]';
b3 = [0, 0, -1]';
cross(b1,b2)==b3

% b/c/d [ab ai bi] to-from
dcmAI = TP.DCM([a1 a2 a3], "A", "I");
dcmBI = TP.DCM([b1 b2 b3], "B", "I");

dcmBA = ROT.DCMmul(dcmBI, ROT.DCMtrs(dcmAI));
% g
dcmAB = ROT.DCMtrs(dcmBA);
% f
dcmAA = ROT.DCMmul(dcmAB,dcmBA);
% h NO

% 2 skewsym

% 3 (WIP)


% 4
eu21 = TP.EulerAng(deg2rad([-15, 25, 10]), [3,2,1], "A", "N");
% a
dcm21 = ROT.euler2dcm(eu21);
% b/c
prv21 = ROT.dcm2PRV(dcm21);
% d
dcm21.Mat*prv21.x;
% e
q21 = ROT.PRV2EP(prv21);
% f YES
norm(q21.x);

% 5
euBN = TP.EulerAng(deg2rad([60, -45, 30]), [3,2,1], "B", "N");

dcmAN = ROT.euler2dcm(eu21);
dcmBN = ROT.euler2dcm(euBN);

dcmBA = ROT.DCMmul(dcmBN, ROT.DCMtrs(dcmAN));
ROT.dcm2euler(dcmBA,'XYZ');

% 6
prv6 = TP.PRV(deg2rad(45.), [1, 1, 1]./sqrt(3), "B", "N");
dcm6 = ROT.PRV2dcm(prv6);
ROT.dcm2euler(dcm6,'XYZ');
ROT.dcm2euler(ROT.DCMtrs(dcm6),'XYZ');