include("Newton.jl")

function ea(h::Float64, u::Vector{Float64}, th::Float64, f::Function)   
    return u + f(u, th)*h
end

function ei(h::Float64, u::Vector{Float64}, th::Float64, f::Function)   
    u1::Vector{Float64} = u
    function F(y::Vector{Float64}, t::Float64)
        return (y - u - h*f(y, th + h))
    end
    u1 = newton(u1, F)
    return u1
end

function cn(h::Float64, u::Vector{Float64}, th::Float64, f::Function)   
    u1::Vector{Float64} = u
    function F(y::Vector{Float64}, t::Float64)
        return (y - u - h/2*(f(y, th + h) + f(u, th)))
    end
    u1 = newton(u1, F)
    return u1
end

function hn(h::Float64, u::Vector{Float64}, th::Float64, f::Function)   
    uu::Vector{Float64} = u + f(u, th)*h
    return u + h/2*(f(u, th) + f(uu, th))
end

function rk4(h::Float64, u::Vector{Float64}, th::Float64, f::Function)   
    k::Matrix{Float64} = zeros(length(u), 4)
    k[:,1] = f(u, th)
    k[:,2] = f(u + k[:,1]*h/2, th + h/2)
    k[:,3] = f(u + k[:,2]*h/2, th + h/2)
    k[:,4] = f(u + k[:,3]*h, th + h)
    return u + h/6*(k*[1,2,2,1])
end

# function am(h, u, th, f, us1, us2)
#     uu = u + f(u, th)*h*23/12 - f(us1, th - h)*h*4/3 + f(us2, th - 2*h)*h*5/12
#     return u + f(u, th)*h*19/24 - f(us1, th - h)*h*5/24 + f(us2, th - 2*h)*h*1/24 + f(uu, th + h)*h*3/8
# end

# function amx(h, u, th, f, us)
#     uu = u + f(u, th)*h*1901/720 - f(us[0], th - h)*h*1387/360 + f(us[1], th - 2*h)*h*109/30 - f(us[2], th - 3*h)*h*637/360 + f(us[3], th - 4*h)*h*251/720
#     #uu=u+h*[1901/720, 1387/360, 109/30, 637/360, 251/720]@[f(u,th), f(us[0],th-h), f(us[1],th-2*h), f(us[2],th-3*h), f(us[3],th-4*h)]
#     return u + f(u, th)*h*646/720 - f(us[0], th - h)*h*264/720 + f(us[1], th - 2*h)*h*106/720 - f(us[2], th - 3*h)*h*19/720 + f(uu, th + h)*h*251/720
# end

# function lf(h, u, th, f, us1, us2)
#     return us1 + f(u, th)*h*2
# end