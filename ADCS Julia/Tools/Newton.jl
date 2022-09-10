function newton(u1::Float64, F::Function; max_iter::Int=100, eps::Float64=1e-14, eps_jac::Float64=1e-8)
    p::Int64 = 0
    while norm(F(u1)) > eps && p < max_iter
        p = p + 1
        A::Float64 = ( F(u1+eps_jac) - F(u1-eps_jac) )/(2*eps_jac)
        u1 = u1 - F(u1)/A
    end
    println("Newton iterations: ", p)
    return u1
end


function newton(u1::Vector{Float64}, F::Function; max_iter::Int=100, eps::Float64=1e-14, eps_jac::Float64=1e-8)
    p::Int64 = 0
    while norm(F(u1)) > eps && p < max_iter
        p = p + 1
        u1 = u1 - Jc(F, u1, eps=eps_jac)\F(u1)
    end
    println("Newton iterations: ", p)
    return u1
end

function Jc(F::Function, U0::Vector{Float64}; eps::Float64=1e-8)
    N = length(U0)
    A = Matrix{Float64}(undef, N, N)
    delta = Vector{Float64}(undef, N)
    for j = 1:N
            delta[:] .= 0
            delta[j] = eps
            A[:,j] = ( F(U0+delta) - F(U0-delta) )/(2*eps)
    end
    return A
end