using Plots
k = 1.380649e-23
T = 273
m = 1.6735e-27

function f(v,T)
    return sqrt.((m/(2*pi*k*T))^3)*4*pi.*v.^2 .*exp.(-m.*v.^2 ./(2*k*T))
end

v = 1:10000
using Plots
plot(f(v,T))

for T = 200:100:1000
    plot!(f(v,T))
end
display(plot!())