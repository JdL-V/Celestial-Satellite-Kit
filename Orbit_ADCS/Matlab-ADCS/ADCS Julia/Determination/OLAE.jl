include("../types.jl")

function OLAE(vkb::Matrix{Float64}, vkn::Matrix{Float64}, weight::Vector{Float64})
    d::Vector{Float64} = vec(vkb - vkn)
    N::Int64 = length(d)
    S::Matrix{Float64} = zeros(N,3)
    for i = 1:Int(N/3)
        S[i:i+2,:] = cross2mat(vkb[:,i] + vkn[:,i])
    end
    W::Matrix{Float64} = Matrix(I,N,N)
    qcr = CRP(inv(S'*W*S)*S'*W*d, "B", "N")

    return CRP2dcm(qcr)
end