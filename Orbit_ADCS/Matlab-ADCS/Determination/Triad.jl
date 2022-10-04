function triad(v1b::Vector{Float64}, v2b::Vector{Float64}, v1n::Vector{Float64}, v2n::Vector{Float64})
    bt::Matrix{Float64} = triadbasegen(v1b, v2b)
    nt::Matrix{Float64} = triadbasegen(v1n, v2n)
    return DCM(bt*nt', "B", "N")
end

function triadbasegen(v1::Vector{Float64}, v2::Vector{Float64})
    t::Matrix{Float64} = zeros(3,3)
    t[:,1] = v1
    t[:,2] = normalize(cross(v1, v2))
    t[:,3] = cross(t[:,1], t[:,2])
    return t
end