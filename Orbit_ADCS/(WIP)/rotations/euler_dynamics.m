function du = euler_dynamics(u, M, I0)
% rotations
i = [1 0 0]';
j = [0 1 0]';
k = [0 0 1]';

I = u(1:3);
J = u(4:6);
K = u(7:9);
w = u(10:12);

A = I0(1,1);
B = I0(2,2);
C = I0(3,3);

dI = cross(w, I);
dJ = cross(w, J);
dK = cross(w, K);
    
dw = -[(C - B)/A*w(2)*w(3);
          (A - C)/B*w(1)*w(3);
          (B - A)/C*w(2)*w(1)];

du = [dI; dJ; dK; dw + M./[A B C]'];
end

