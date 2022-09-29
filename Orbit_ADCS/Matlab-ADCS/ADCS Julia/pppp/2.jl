include("../AstroDynamics/astroConstants.jl")
Re, mu , J2= astroConstants([23, 13, 9])
rp = Re + 300
ra = Re + 10000

a = (rp + ra)/2
e = (ra - rp)/(ra + rp)
p = a*(1 - e^2)
vp = sqrt(mu/p)*(1 + e)
va = sqrt(mu/p)*(1 - e)
T = 2*pi*sqrt(a^3/mu)/3600

display([a,e,vp,va,T])

include("../Tools/Newton.jl")
rp = Re + 800
v = 8

function FFF(x)
    return [v - sqrt(mu/x[1])*(1 + x[2]), rp - x[1]/(1+ x[2])]
end

p, e = newton([rp, 0.5], FFF)
ra = p/(1 - e)
a = (rp + ra)/2
hmax = ra - Re

display([a, e, ra, hmax])
OM = 2*pi/(365.25*24*3600)

a = [1, 1.5, 2]*1e4

ee(a) = sqrt.(1 .- sqrt.(3*J2*Re^2*sqrt.(mu./(5*a.^7))/(2*OM)))
aa(a) = 1 - sqrt.(3*J2*Re^2*sqrt(mu/(5*a^7))/(2*OM))
rpp(a) = a.*(1 .- ee(a))

ee(1e4)
ee(1.5e4)
ee(2e4)

ee(a)

a0 = newton(1e4, aa)
display(a0)
rpp(a0 + eps(a0)*1e2)
a = (a0 + 1e-10):10:1e5
using Plots
display(plot(rpp(a), a))

rpp2(a) = a.*(1 .- ee(a)) - 200-Re
rp0 = newton(1e4, rpp2, max_iter=1000, eps_N=1e-4, eps_jac=1e-4)