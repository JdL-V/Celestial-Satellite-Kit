include("../types.jl")
include("../RefRotations.jl")
include("../RefDetermination.jl")
include("../Tools/Cauchy.jl")
include("../Tools/Newton.jl")

M = 1e3
At = 1e4
ps = 2*1361/299792458
as = 2*ps*At/M
f(t) = as/2*t^2 - 4e8
display(newton(1e6, f, max_iter=1000, eps_N=1e-6, eps_jac=1e-8))

# function fF5(u,p,t)
#     du = [u[2], as]
#     return du
# end
# u0 = [0,0]
# @time begin
#     using DifferentialEquations
#     condition(u,t,integrator) = u[1] > 384.4e6
#     function dose!(integrator)     
#         terminate!(integrator)
#     end
#     cb = cb = DiscreteCallback(condition, dose!, save_positions=(false,false))
#     a = ODEProblem(fF5, u0, (0.,1e15), callback=cb)
#     s1 = solve(a, DP8(), reltol=1e-14, abstol=1e-14)
# end

display(plot(s1, layout=(2,1)))

####################################

I3 = rand()*100
I2 = I1 = rand()*100
b = rand()/5

Ic = [I1 0 0; 0 I2 0; 0 0 I3]
function fF3(u,p,t)
    Lc = -b.*u
    du = -cross2mat(u)*Ib*u+ Lc
    return du
end
u0 = [rand(), rand(), rand()]

@time begin
    using DifferentialEquations
    a = ODEProblem(fF3, u0, (0.,100.))
    s = solve(a, DP8(), reltol=1e-14, abstol=1e-14)
end

display(plot(s))

# dw3 = -b*w3/I33
# w3 = K*e^(-b*t/I33)