# 0.45 86k

@time begin
include("RefRotations.jl")
include("Tools/Cauchy.jl")
include("Tools/MathTools.jl")
include("RefDetermination.jl")
include("types.jl")

# using NLsolve
# using ReferenceFrameRotations
# using DifferentialEquations
# using LinearAlgebra

eu1 = EulerAng([90., 90.], [3, 1], "B", "N")
eu2 = EulerAng([45., -45.], [3, 2], "B", "N")

A = euler2dcm(eu1).Mat*euler2dcm(eu2).Mat
eu3 = EulerAng([45., -30., 0.], [3, 2, 3], "B", "N")
B = euler2dcm(eu3).Mat[[2,3,1],:]

eu4 = EulerAng([20., 10., -10.], [3, 2, 1], "B", "N")
C = euler2dcm(eu4)
display(C)

s = dcm2euler(C, :ZXZ)
eus = EulerAng(rad2deg.(s.Angles), [3, 1, 3], "B", "N")
C = euler2dcm(eus);
display(C)
display(rad2deg.(s.Angles))

eu5 = EulerAng([10., 20., 30.], [3, 2, 1], "B", "B")
D = euler2dcm(eu5)
display(rad2deg.(dcm2euler(D, :ZXZ).Angles))

eu6 = EulerAng([-5., 5., 5.], [3, 2, 1], "B", "B")
E = dcmmul(D, euler2dcm(eu6)) # DCM(D.Mat*(euler2dcm(eu6).Mat)', "B", "N")
display(rad2deg.(dcm2euler(E, :ZYX).Angles))

F = dcmmul(E, euler2dcm(eu6)).Mat - D.Mat
display(F)

function ff0(u,t)
    # 321
    w(t) = deg2rad(20).*[sin(0.1*t), 0.01, cos(0.1*t)]
    M::Matrix{Float64} = om2euler([u[1], u[2], u[3]], :ZYX)
    du::Vector{Float64} = M*w(t)
    return du
end
u0 = deg2rad.([40, 30, 80])
tt, uu = ca(Int(1e5), 42., u0, ff0, rk4)
norm(uu[end,:])

G = [0.925417   0.0296956    0.377786
0.336824   -0.521281    -0.784102  
0.173648    0.852869    -0.492404]'

eug = EulerAng([20., -10., 120.], [3, 2, 1], "B", "N")
G = euler2dcm(eug)

G = [1 0 0; 0 0 1; 0 -1 0]
G = G
display(G)
    
ϕ = acos(0.5*(G[1,1]+G[2,2]+G[3,3]-1))
display(ϕ)
e = 1/(2*sin(ϕ))*[G[2,3]-G[3,2] G[3,1]-G[1,3] G[1,2]-G[2,1]]

prv = PRV(pi/4, [1., 0., 0.], "B", "N")

A = PRV2om(prv)
m = A*[3., 1., 2.]
B = om2PRV(prv)
display(B*m)
q = quaternion([0.235702, 0.471405, -0.471405, 0.707107], "B", "N")
b = EP2dcm(q)

display(dcm2EP(b))

a = dcm2EP(DCM(trs([-0.529403 -0.474115  0.703525
        -0.467056 -0.529403 -0.708231 
         0.708231 -0.703525  0.0588291]), "B", "N"))
display(a)
display(sum(a.x.^2))

function ff1(u,t)
    # 321
    w(t) = deg2rad(20).*[sin(0.1*t), 0.01, cos(0.1*t)]
    M::Matrix{Float64} = om2EP(u)
    du::Vector{Float64} = M*w(t)
    return du
end
u0 = [0.408248, 0., 0.408248, 0.816497]
tt, uu = ca(Int(1e5), 42., u0, ff1, rk4)
norm(uu[end,2:4])

m = CRP([0.1, 0.2, 0.3], "B", "N")

display(CRP2dcm(m))

A = DCM([0.333333    -0.666667    0.666667
0.871795    0.487179    0.0512821
-0.358974    0.564103    0.74359], "B", "N")
h = dcm2CRP(A)
display(h)
h.x = -h.x
display(CRP2dcm(h))

function ff2(u,t)
    # 321
    w(t) = deg2rad(3).*[sin(0.1*t), 0.01, cos(0.1*t)]
    M::Matrix{Float64} = om2CRP(u)
    du::Vector{Float64} = M*w(t)
    return du
end
u0 = [0.4, 0.2, -0.1]
tt, uu = ca(Int(1e5), 42., u0, ff2, rk4)
norm(uu[end,:])

    A = DCM([0.925417   0.0296956    0.377786
    0.336824   -0.521281    -0.784102  
    0.173648    0.852869    -0.492404],"B","N")
display(dcm2EP(A))
display(sheppard(A))

mrp1 = MRP([0.1, 0.2, 0.3], "B", "N")
mrp2 = MRP([-0.1, 0.3, 0.1], "B", "N")
mrp3 = MRP(-[0.5, 0.3, 0.1], "B", "N")
MRP2dcm(mrp1)

A = DCM([0.763314    0.0946746    -0.639053
-0.568047    -0.372781    -0.733728
-0.307692    0.923077    -0.230769],"B","N")
dcm2MRP(A)

display(SumMRP(mrp1, mrp2))
display(SumMRP(mrp3, mrp1))


function ff3(u,t)
    321
    w(t)::Vector{Float64} = deg2rad(20).*[sin(0.1*t), 0.01, cos(0.1*t)]
    M::Matrix{Float64} = om2MRP(u)
    du::Vector{Float64} = M*w(t)
    return du
end
u0 = [0.4, 0.2, -0.1]
tt, uu = MRPdiff(Int(1e5), 42., u0, ff3, rk4)
norm(uu[end,:])


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

AW::Matrix{Float64} = [    0.969846    0.17101    0.173648
-0.200706    0.96461    0.17101
-0.138258    -0.200706    0.969846]
BW::Matrix{Float64} = [    0.963592    0.187303    0.190809
-0.223042    0.956645    0.187303
-0.147454    -0.223042    0.963592]

ϕW::Float64 = rad2deg(dcm2PRV(DCM(AW*BW',"B","N")).x[1])

end