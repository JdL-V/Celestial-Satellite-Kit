%1
% I1 > I3
% m/12*H^2 + m/4*R^2 > m/2*R^2
% H = sqrt(3)*R

%2
Ib = [100   0       0;
      0     120     0;
      0     0       80];
RoG = [0 0 7e3]';
mu = astroConstants(13);
euBG = EulerAng(ones(1,3).*pi/4, [3,2,1], "B", "G");
dcmBG = euler2dcm(euBG);
RoB = dcmBG.Mat*RoG;

T = 3*mu/norm(RoB)^5*skewsym(RoB)*Ib*RoB

%3
I3 = [8     0       0;
      0     12      0;
      0     0       10];
In = diag(I3);
wz = 0.1;
cond1 = (In(2) - In(3))*wz
cond2 = (In(1) - In(3))*wz

%4
%Dual spin stability schaub