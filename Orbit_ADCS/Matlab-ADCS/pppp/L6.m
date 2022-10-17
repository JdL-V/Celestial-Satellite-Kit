ROT = RotModule;
TP = types;
format compact

a1 = [0.5, sqrt(3)/2, 0]';
a2 = [0., 0., 1.]';
a3 = [sqrt(3)/2, -0.5, 0]';
cross(a1,a2)==a3

b1 = [0, 1, 0]';
b2 = [1, 0, 0]';
b3 = [0, 0, -1]';
cross(b1,b2)==b3

% to-from
dcmAI = TP.DCM([a1 a2 a3], "A", "I")
dcmBI = TP.DCM([b1 b2 b3], "B", "I")

dcmBA = ROT.DCMmul(dcmBI, ROT.DCMtrs(dcmAI))
dcmAB = ROT.DCMtrs(dcmBA)
dcmAA = ROT.DCMmul(dcmAB,dcmBA)


eu21 = TP.EulerAng(deg2rad([-15, 25, 10]), [3,2,1], "A", "N");
dcm21 = ROT.euler2dcm(eu21);
prv21 = ROT.dcm2PRV(dcm21);
dcm21.Mat*prv21.x;
q21 = ROT.PRV2EP(prv21);

norm(q21.x);

euBN = TP.EulerAng(deg2rad([60, -45, 30]), [3,2,1], "B", "N");

dcmAN = ROT.euler2dcm(eu21);
dcmBN = ROT.euler2dcm(euBN);

dcmBA = ROT.DCMmul(dcmBN, ROT.DCMtrs(dcmAN));
ROT.dcm2euler(dcmBA,'XYZ');

prv6 = TP.PRV(deg2rad(45.), [1, 1, 1]./sqrt(3), "B", "N");
dcm6 = ROT.PRV2dcm(prv6);
ROT.dcm2euler(dcm6,'XYZ');
ROT.dcm2euler(ROT.DCMtrs(dcm6),'XYZ');