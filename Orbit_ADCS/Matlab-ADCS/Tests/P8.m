IB1 = [1500. 0. -1000.; 0. 2700. 0.; -1000. 0. 3000.];
[I01, dcm0B1] = getPrincipalInertia(IB1)

IB = [20. -10. -15.; -10. 30. 0.; -15. 0. 15.];
[I0, dcm0B] = getPrincipalInertia(IB)

T = diag(I0)'

getRotKineticE(IB, [1., 2., 3.])
getRotKineticE(I0, dcm0B*[1., 2., 3.]')

% 0-2π/dθ  *  0-R/kr^3*r^4 dr  *  0-π/sin^3(ϕ) dϕ pag 478 Bramanti AG2
% kR^5/4 

% Ixx = /y^2 + z^2 dm
% Iyy = /x^2 + z^2 dm
% Izz = /y^2 + x^2 dm
% Ixx + Izz - Iyy = /2y^2 dm