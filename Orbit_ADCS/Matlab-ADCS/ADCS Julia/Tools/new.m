function [u1] = new(u1, F)
p = 0;
while norm(F(u1)) > 1e-15 && p < 50
    p = p + 1;
    u1 = u1 - Jc(F, u1, 1e-8)\F(u1);
end
end

function [A] = Jc(F, U0, eps)
N = size(U0,1);
A = zeros(N,N);
delta = zeros(N,1);
for j=1:N
    delta(:) = 0;
    delta(j) = eps;
    A(:,j) = ( F(U0+delta) - F(U0-delta) )/(2*eps);
end
end