function var = eznewton(u1, F, max_iter)
    p = 0;
    if nargin == 2
        max_iter = 100;
    end

    while norm(F(u1)) > 1e-14 && p < max_iter
        p = p + 1;
        u1 = u1 - Jc(F, u1, eps(u1)*1e8)\F(u1);
    end
    fprintf("Newton iterations: %d\n", p)
    var = u1;
end

function var = Jc(F, U0, eps_J)
    N = length(U0);
    A = zeros(N, N);
    delta = zeros(1, N);
    for j = 1:N
        delta(:) = 0;
        delta(j) = eps_J;
        A(:,j) = ( F(U0+delta) - F(U0-delta) )/(2*eps_J);
    end
    var =  A;
end