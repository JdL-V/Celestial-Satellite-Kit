@time begin

include("RefDynamics.jl")
include("RefRotations.jl")
include("Tools/MathTools.jl")
include("Tools/Cauchy.jl")
include("types.jl")

M = [1., 1., 2., 2.]
R = [[1., -1., 2.], [-1., -3., 2.], [2., -1., -1.], [3., -1., -2.]]
dR = [[2., 1., 1.], [0., -1., 1.], [3., 2., -1.], [0., 0., 1.]]

display(getCM(M,R))
display(getKineticE(M, dR))
display(getAngularMom(M,R,dR,[0.,0.,0.],[0.,0.,0.]))
display(getAngularMom(M,R,dR,getCM(M,R),getCM(M,dR)))

eu = EulerAng(deg2rad.([-10.,10.,5.]), [3,2,1], "B", "N")
dcm = euler2dcm(eu)
Ib = [10. 1. -1.; 1. 5. 1.; -1. 1. 8.]

Ib*(dcm.Mat*[0.01,-0.01,0.01])

M = [12.5]
In = rotateInertiaFrame(Ib, trs(dcm.Mat))
In = translateInertiaFrame(In, M, [-0.5,0.5,0.25])

Ibp = translateInertiaFrame(Ib, M, dcm.Mat*[-0.5,0.5,0.25])

Ibb = getPrincipalInertia(Ib)
display(Ibb[1])
display(Ibb[2])

mrp = MRP([0.1, 0.2, 0.3], "B", "N")
In = rotateInertiaFrame(Ib, MRP2dcm(mrp).Mat)

display(getRotKineticE(Ib,[0.01, -0.01, 0.01]))

function fF(u,t)
    Lc = [0,1,1]
    du = [om2MRP(u[1:3])*u[4:6]; -cross2mat(u[4:6])*Ib*u[4:6] + Lc]
    return du
end
u0 = [0.4, 0.2, -0.1, 3, 4, 5]
tt, uu = MRPdiff(Int(1e5), 42., u0, fF, rk4)
norm(uu[end,:])


end