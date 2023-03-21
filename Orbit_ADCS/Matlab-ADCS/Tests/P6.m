format compact

% 1
% a Check orthonormal reference frame
a1 = [0.5, sqrt(3)/2, 0]';
a2 = [0., 0., 1.]';
a3 = [sqrt(3)/2, -0.5, 0]';
(cross(a1,a2)==a3)'

b1 = [0, 1, 0]';
b2 = [1, 0, 0]';
b3 = [0, 0, -1]';
(cross(b1,b2)==b3)'

% b/c/d [ab ai bi] to-from
dcmAI = DCM([a1 a2 a3]', "A", "I");
dcmBI = DCM([b1 b2 b3]', "B", "I");

dcmBA = DCMmul(dcmBI, DCMtrs(dcmAI));
% g
dcmAB = DCMtrs(dcmBA);
% f
dcmAA = DCMmul(dcmAB,dcmBA);
% h NO

% 4
eu21 = EulerAng(deg2rad([10, 25, -15]), [3,2,1], "A", "N");
% a
dcm21 = euler2dcm(eu21);
% b/c
prv21 = dcm2PRV(dcm21);
% d
dcm21.Mat*prv21.x;
% e
q21 = PRV2EP(prv21);
% f YES
norm(q21.x);
tic
DDD = dcm2EP(dcm21);
toc
tic
DDD2 = sheppard(dcm21);
toc
% 5
euAN = EulerAng(deg2rad([30, -45, 60]), [3,2,1], "B", "N");

dcmBN = euler2dcm(eu21);
dcmAN = euler2dcm(euAN);

dcmAB2 = DCMmul(dcmAN, DCMtrs(dcmBN));
eu5 = dcm2euler(dcmAB2,'ZYX');

% 6
prv6 = PRV(deg2rad(45.), [1, 1, 1]./sqrt(3), "B", "N");
dcm6 = PRV2dcm(prv6);
eu61 = dcm2euler(dcm6,'ZYX');
eu62 = dcm2euler(DCMtrs(dcm6),'ZYX');