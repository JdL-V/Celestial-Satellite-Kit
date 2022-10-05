function var = newton(u1, F, max_iter, eps_N, eps_jac)
    p = 0;
    while norm(F(u1)) > eps_N && p < max_iter
        p = p + 1;
        u1 = u1 - Jc(F, u1, eps_jac)\F(u1);
    end
    fprintf("Newton iterations: %d", p)
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