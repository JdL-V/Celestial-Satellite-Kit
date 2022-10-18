prvBN = PRV(rad2deg(45), [0.,0.,1.], "B", "N")
dcmBN = PRV2dcm(prvBN)
neN = [-1.,0.,0.]
nsN = [0.,1.,0.]
neB = dcmBN.Mat*neN
nsB = dcmBN.Mat*nsN
dcmBN2 = triad(neB, nsB, neN, nsN)
display(dcmBN.Mat)
display(dcmBN2.Mat)
dcmBN3 = qMethod([neB nsB], [neN nsN], [1.,1.])
dcmBN4 = QUEST([neB nsB], [neN nsN], [1.,1.])
display(dcmBN3.Mat)
display(dcmBN4.Mat)


dcmBN = triad([0.8273, 0.5541, 0.0920], [-0.8285, 0.5522, -0.0955], [-0.1517, -0.9669, 0.2050], [-0.8393, 0.4494, -0.3044])

function feuler(u,t)
    # 321
    w(t) = deg2rad(5).*[sin(0.1*t), 0.01, cos(0.1*t)]
    M::Matrix{Float64} = om2euler([u[1], u[2], u[3]], :ZYX)
    du::Vector{Float64} = M*w(t)
    return du
end
u0 = deg2rad.([80, 30, 40])
tt, uu = ca(Int(1e5), 60., u0, feuler, rk4)
uu[end,:]


# function fq(u,t)
#     # 321
#     w(t) = deg2rad(50).*[sin(0.1*t), 0.01, cos(0.1*t)]
#     M::Matrix{Float64} = om2EP(u)
#     du::Vector{Float64} = M*w(t)
#     return du
# end

# u0 = dcm2EP(euler2dcm(EulerAng(u0,[3,2,1],"B","N"))).x
# @time begin
# tt, uu = EPdiff(Int(1e5), 60., u0, fq, rk4)
# end
# display(uu[end,:])

function fq2(u,p,t)
    # 321
    w(t) = deg2rad(50).*[sin(0.1*t), 0.01, cos(0.1*t)]
    M::Matrix{Float64} = om2EP(u)
    du::Vector{Float64} = M*w(t)
    return du
end

@time begin
u0 = [1, 0, 0, 0]
using DifferentialEquations
condition(u,t,integrator) = norm(u) != 1
function dose!(integrator)     
   integrator.u = normalize(integrator.u)
end
cb = cb = DiscreteCallback(condition, dose!, save_positions=(false,true))
a = ODEProblem(fq2, u0, (0.,60.))
s = solve(a, DP8(), reltol=1e-14, abstol=1e-14, callback=cb)
display(norm.(s.u))
display(minimum(norm.(s.u)))
end


using Plots

display(plot(tt, uu))
display(plot(s))