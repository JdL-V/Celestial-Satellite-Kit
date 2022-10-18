include("../RefDynamics.jl")
include("../RefRotations.jl")
include("../Tools/MathTools.jl")
include("../Tools/Cauchy.jl")
include("../types.jl")

IB = [1500. 0. -1000.; 0. 2700. 0.; -1000. 0. 3000.]
I0, dcm0B = getPrincipalInertia(IB)

IB = [20. -10. -15.; -10. 30. 0.; -15. 0. 15.]
I0, dcm0B = getPrincipalInertia(IB)

T = [I0[1,1], I0[2,2], I0[3,3]]

display(getRotKineticE(IB, [1., 2., 3.]))
display(getRotKineticE(I0, dcm0B*[1., 2., 3.]))

# 0-2π/dθ  *  0-R/kr^3*r^4 dr  *  0-π/sin^3(ϕ) dϕ pag 478 Bramanti AG2
# kR^5/4 

# Ixx = /y^2 + z^2 dm
# Iyy = /x^2 + z^2 dm
# Izz = /y^2 + x^2 dm
# Ixx + Izz - Iyy = /2y^2 dm