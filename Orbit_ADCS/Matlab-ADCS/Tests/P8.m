%1
IB1 = [1500. 0. -1000.; 0. 2700. 0.; -1000. 0. 3000.];
[I01, dcm0B1] = getPrincipalInertia(IB1)

%2
IB = [20. -10. -15.;
     -10.  30.  0.; 
     -15.  0.   15.];
[I0, dcm0B] = getPrincipalInertia(IB)
dcm0B*IB*dcm0B' 

T = diag(I0)'/2

getRotKineticE(IB, [1., 2., 3.])
getRotKineticE(I0, dcm0B*[1., 2., 3.]')

%3
% 0-2π// dθ  *  0-R// rho(kr^3)*r^4 dr  *  0-π// sin^3(ϕ) dϕ pag 478
% kπ/3*R^8

%4
% Ixx = /y^2 + z^2 dm
% Iyy = /x^2 + z^2 dm
% Izz = /y^2 + x^2 dm
% Ixx + Izz - Iyy = /2y^2 dm















% Ib = [150. 0. -100.; 0. 250. 0.; -100. 0. 300.];
% [I0a, dcm0Ba] = getPrincipalInertia(Ib)
% dcm0Ba*Ib*dcm0Ba' 

% neB = normalize([0 1  -1],'norm')'
% nsB = normalize([0 1 1],'norm')'

% neN = [0 0 1]'
% nsN = [0 1 0]'

% dcmBN2 = qMethod([neB nsB], [neN nsN], [1,1])