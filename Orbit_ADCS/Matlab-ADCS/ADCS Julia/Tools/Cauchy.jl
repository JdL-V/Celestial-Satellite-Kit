include("ODE.jl")

function ca(nn::Int64, tmax::Float64, u0::Vector{Float64}, f::Function, sch::Function)
    h::Float64  = tmax/(nn-1)
    th::LinRange{Float64, Int64} = LinRange(0, tmax, nn)
    u::Matrix{Float64} = zeros(nn, length(u0))
    u[1,:] = u0

    for n::Int64 = 1:(nn-1)
        u[n+1,:] = sch(h, u[n,:], th[n], f)
    end
    return th, u
end