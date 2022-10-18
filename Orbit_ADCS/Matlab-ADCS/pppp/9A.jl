include("../RefDynamics.jl")
include("../RefRotations.jl")
include("../Tools/MathTools.jl")
include("../Tools/Cauchy.jl")
include("../types.jl")

using DifferentialEquations
using Plots

m = 1.
R = 1.
Om0 = 0.5
Ib = [4*m*R^2/3 + m*R^2/4 0 0; 0 4*m*R^2/3 + m*R^2/4 0; 0 0 m*R^2/2]

# function fF(u,t)
#     Lc = [0,0,0]
#     du = [om2MRP(u[1:3])*u[4:6]; -skewsym(u[4:6])*Ib*u[4:6] + Lc]
#     return du
# end
# u0 = [0., 0., 0., 3/5*Om0, 0., 4/5*Om0]
# tt, uu = MRPdiff(Int(1e5), 60., u0, fF, rk4)
# uu[end,:]

# using Plots
# ymax = maximum(uu)
# ymin = minimum(uu)
# dy = ymax - ymin
# ymax = ymax + dy*0.1
# ymin = ymin - dy*0.1
# display(plot(tt,uu[:,4:6], ylim=(ymin, ymax)))
# display(plot(tt,uu[:,1:3], ylim=(ymin, ymax)))

#########################################

function fF3(u,p,t)
    Lc = [0,0,0]
    du = [om2MRP(u[1:3])*u[4:6]; -skewsym(u[4:6])*Ib*u[4:6] + Lc]
    return du
end
u0 = [0., 0., 0., 3/5*Om0, 0., 4/5*Om0]

@time begin
    using DifferentialEquations
    condition(u,t,integrator) = dot(u[1:3], u[1:3]) > 1
    function dose!(integrator)     
       integrator.u[1:3] = -integrator.u[1:3]./dot(integrator.u[1:3],integrator.u[1:3])
    end
    cb = cb = DiscreteCallback(condition, dose!, save_positions=(true,true))
    a = ODEProblem(fF3, u0, (0.,60.))
    s = solve(a, DP8(), reltol=1e-14, abstol=1e-14, callback=cb)
end

# function f2F(u,p,t)
#     du = -skewsym(u)*Ib*u
#     return du
# end

# a = ODEProblem(f2F, u0, (0.,60.))
# s = solve(a, DP8(), reltol=1e-8, abstol=1e-8)

display(plot(s, vars=[4,5,6]))
display(plot(s, vars=[1,2,3]))

#####################################


I2 = [98. 0. 0.; 0. 102. 0.; 0. 0. 150.]
w = [0.1, 0.02, 0.5]
H = I2*w
β = acosd(dot(normalize(H),normalize(w)))
nut = acosd(dot(normalize(w),[1,0,0]))

function fF4(u,p,t)
    Lc = [0,0,0]
    du = [om2EP(u[1:4])*u[5:7]; -skewsym(u[5:7])*Ib*u[5:7] + Lc]
    return du
end
u0 = [1., 0., 0., 0., w[1], w[2], w[3]]

@time begin
    using DifferentialEquations
    condition(u,t,integrator) = norm(u) != 1
    function dose!(integrator)     
       integrator.u[1:4] = normalize(integrator.u[1:4])
    end
    cb = cb = DiscreteCallback(condition, dose!, save_positions=(true,true))
    a = ODEProblem(fF4, u0, (0.,60.))
    s = solve(a, DP8(), reltol=1e-14, abstol=1e-14, callback=cb)
    display(norm.(s.u))
    display(minimum(norm.(s.u)))
end

display(plot(s, vars=[5,6,7]))
display(plot(s, vars=[1,2,3,4]))

for i = 1:111:553
w = s.u[i][5:7]
H = I2*w
β = acosd(dot(normalize(H),normalize(w)))
nut = acosd(dot(normalize(w),[0,0,1]))
display(β)
display(nut)
end


#################################################

I2 = [210. 0. 0.; 0. 200. 0.; 0. 0. 118.]
w = [0.05, 0.02, -0.02]

function fF5(u,p,t)
    Lc = [0,0,0]
    du = [om2euler(u[1:3], :ZYX)*u[4:6]; -skewsym(u[4:6])*Ib*u[4:6] + Lc]
    return du
end
u0 = [0., 0., 0., w[1], w[2], w[3]]

@time begin
    using DifferentialEquations
    a = ODEProblem(fF5, u0, (0.,100.))
    s = solve(a, DP8(), reltol=1e-14, abstol=1e-14)
    display(norm.(s.u))
    display(minimum(norm.(s.u)))
end

display(plot(s, vars=[4,5,6]))
display(plot(s, vars=[1,2,3]))