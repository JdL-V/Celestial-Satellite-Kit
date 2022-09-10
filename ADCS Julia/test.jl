@time begin
include("RefRotations.jl")
include("Tools/Cauchy.jl")
include("RefDetermination.jl")

# using NLsolve
# using ReferenceFrameRotations
# using DifferentialEquations
# using LinearAlgebra


# A = euler2dcm(3, 1, 3)(90., 90., 0.)*euler2dcm(3, 2, 3)(45., -45., 0.)
# B = euler2dcm(3, 2, 3)(45., -30., 0.)[[2,3,1],:]

# C = euler2dcm(3, 2, 1)(20., 10., -10.)
# display(C)

# s = dcm2euler(C, :ZXZ)
# C = euler2dcm(3, 1, 3)(rad2deg(s[1]), rad2deg(s[2]), rad2deg(s[3]));
# display(C)
# display(rad2deg.(s))

# D = euler2dcm(3, 2, 1)(10., 20., 30.)
# display(rad2deg.(dcm2euler(D, :ZXZ)))

# E = euler2dcm(3, 2, 1)(10., 20., 30.)*euler2dcm(3, 2, 1)(-5., 5., 5.)'
# display(rad2deg.(dcm2euler(E, :ZYX)))

# F = E*euler2dcm(3, 2, 1)(-5., 5., 5.)-D
# display(F)

# function ff0(u,p,t)
#     # 321
#     w(t) = deg2rad(20).*[sin(0.1*t), 0.01, cos(0.1*t)]
#     M::Matrix{Float64} = om2euler(u[1], u[2], u[3], :ZYX)
#     du::Vector{Float64} = M*w(t)
#     return du
# end
# u0 = deg2rad.([40, 30, 80])
# tt, uu = ca(Int(1e5), 42., u0, ff0, rk4)
# norm(uu[end,:])

# G = [0.925417   0.0296956    0.377786
# 0.336824   -0.521281    -0.784102  
# 0.173648    0.852869    -0.492404]'

# G = euler2dcm(3, 2, 1)(20., -10., 120.)

# G = [1 0 0; 0 0 1; 0 -1 0]
# G = G
# display(G)
    
# ϕ = acos(0.5*(G[1,1]+G[2,2]+G[3,3]-1))
# display(ϕ)
# e = 1/(2*sin(ϕ))*[G[2,3]-G[3,2] G[3,1]-G[1,3] G[1,2]-G[2,1]]

# A = PRV2om(pi/4, [1 0 0])
# m = A*[3 1 2]'
# B = om2PRV(pi/4, [1 0 0])
# display(B*m)

# b = EP2dcm([0.235702 0.471405 -0.471405 0.707107])

# display(dcm2EP(b))

# a = dcm2EP([-0.529403 -0.474115  0.703525
#         -0.467056 -0.529403 -0.708231 
#          0.708231 -0.703525  0.0588291]')
# display(a)
# display(sum(a.^2))

# function ff1(u,t)
#     # 321
#     w(t) = deg2rad(20).*[sin(0.1*t), 0.01, cos(0.1*t)]
#     M::Matrix{Float64} = om2EP(u)
#     du::Vector{Float64} = M*w(t)
#     return du
# end
# u0 = [0.408248, 0., 0.408248, 0.816497]
# tt, uu = ca(Int(1e5), 42., u0, ff1, rk4)
# norm(uu[end,2:4])

# m = [0.1, 0.2, 0.3]

# display(CRP2dcm(m))

# A = [0.333333    -0.666667    0.666667
# 0.871795    0.487179    0.0512821
# -0.358974    0.564103    0.74359]
# h = dcm2CRP(A)
# display(h)
# display(CRP2dcm(-h))

# function ff2(u,t)
#     # 321
#     w(t) = deg2rad(3).*[sin(0.1*t), 0.01, cos(0.1*t)]
#     M::Matrix{Float64} = om2CRP(u)
#     du::Vector{Float64} = M*w(t)
#     return du
# end
# u0 = [0.4, 0.2, -0.1]
# tt, uu = ca(Int(1e5), 42., u0, ff2, rk4)
# norm(uu[end,:])

#     A =[0.925417   0.0296956    0.377786
#     0.336824   -0.521281    -0.784102  
#     0.173648    0.852869    -0.492404]
# display(dcm2EP(A))
# display(sheppard(A))

# MRP2dcm([0.1, 0.2, 0.3])

# A = [0.763314    0.0946746    -0.639053
# -0.568047    -0.372781    -0.733728
# -0.307692    0.923077    -0.230769]
# dcm2MRP(A)

# display(SumMRP([0.1, 0.2, 0.3], [-0.1, 0.3, 0.1]))
# display(SumMRP(-[0.5, 0.3, 0.1], [0.1, 0.2, 0.3]))


# function ff3(u,t)
#     # 321
#     w(t)::Vector{Float64} = deg2rad(20).*[sin(0.1*t), 0.01, cos(0.1*t)]
#     M::Matrix{Float64} = om2MRP(u)
#     du::Vector{Float64} = M*w(t)
#     return du
# end
# u0 = [0.4, 0.2, -0.1]
# tt, uu = MRPdiff(Int(1e5), 42., u0, ff3, rk4)
# norm(uu[end,:])


v1b::Vector{Float64} = [0.8273, 0.5541, -0.0920]
v2b::Vector{Float64} = [-0.8285, 0.5522, -0.0955]
v1n::Vector{Float64} = [-0.1517, -0.9669, 0.2050]
v2n::Vector{Float64} = [-0.8393, 0.4494, -0.3044]
@time begin
    display(triad(v1b, v2b, v1n, v2n))
end
@time begin
    display(qMethod([v1b v2b], [v1n v2n],[1., 1.]))
end
@time begin
    display(QUEST([v1b v2b], [v1n v2n],[1., 1.]))
end
@time begin
    display(OLAE([v1b v2b], [v1n v2n],[1., 1.]))
end

# A::Matrix{Float64} = [    0.969846    0.17101    0.173648
# -0.200706    0.96461    0.17101
# -0.138258    -0.200706    0.969846]
# B::Matrix{Float64} = [    0.963592    0.187303    0.190809
# -0.223042    0.956645    0.187303
# -0.147454    -0.223042    0.963592]

# ϕ::Float64 = rad2deg(dcm2PRV(A*B')[1])

end