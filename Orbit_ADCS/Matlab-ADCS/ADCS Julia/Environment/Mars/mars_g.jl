using LinearAlgebra
function mars_g(r::Float64)
    g::Float64 = 3.721*(3389.5e3./r).^2;
    return g
end

function mars_g(r::Vector{Float64})
    g::Float64 = 3.721*(3389.5e3./norm(r)).^2;
    return g
end