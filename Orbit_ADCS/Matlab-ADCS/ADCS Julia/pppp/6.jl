a1 = [0.5, sqrt(3)/2, 0]
a2 = [0., 0., 1.]
a3 = [sqrt(3)/2, -0.5, 0]
display(cross(a1,a2)==a3)

b1 = [0, 1, 0]
b2 = [1, 0, 0]
b3 = [0, 0, -1]
display(cross(b1,b2)==b3)

# to-from
AI = [a1 a2 a3]
BI = [b1 b2 b3]

BA = BI*AI'
AB = BA'
display(AB*BA)

include("../types.jl")
include("../RefRotations.jl")

eu21 = EulerAng(deg2rad.([-15, 25, 10]), [3,2,1], "A", "N")
dcm21 = euler2dcm(eu21)
prv21 = dcm2PRV(dcm21)
dcm21.Mat*prv21.x
q21 = PRV2EP(prv21)

norm(q21.x)

euBN = EulerAng(deg2rad.([60, -45, 30]), [3,2,1], "B", "N")

dcmAN = euler2dcm(eu21)
dcmBN = euler2dcm(euBN)

dcmBA = dcmmul(dcmBN, dcmtrs(dcmAN))
dcm2euler(dcmBA)

prv6 = PRV(deg2rad(45.), [1, 1, 1]./sqrt(3), "B", "N")
dcm6 = PRV2dcm(prv6)
dcm2euler(dcm6)
dcm2euler(dcmtrs(dcm6))