using LinearAlgebra
function earth_g(r::Float64)
    g::Float64 = 9.81*(6371e3./r).^2
    return g
end

function earth_g(r::Vector{Float64})
    g::Float64 = 9.81*(6371e3./norm(r)).^2
    return g
end