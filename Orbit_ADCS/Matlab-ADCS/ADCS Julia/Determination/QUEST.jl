include("../Rotations/EPRot.jl")
include("../Tools/newton.jl")
include("../types.jl")

function QUEST(vkb::Matrix{Float64}, vkn::Matrix{Float64}, weight::Vector{Float64})
    B::Matrix{Float64} = zeros(3,3)
    for i::Int64 = 1:size(vkb, 2)
        B = B + weight[i]*vkb[:,i]*vkn[:,i]'
    end

    S::Matrix{Float64} = B + B'
    Σ::Float64 = tr(B)
    Z::Vector{Float64} = [B[2,3] - B[3,2], B[3,1] - B[1,3], B[1,2] - B[2,1]]

    K::Matrix{Float64} = [Σ Z'
                          Z S-Σ.*Matrix(I,3,3)]
    
    f(s) = det(K - s[1].*Matrix(I,4,4))
    λ0::Float64 = sum(weight)
    λopt::Float64 = newton(λ0, f)

    qcr = CRP(inv((λopt + Σ)*Matrix(I,3,3) - S)*Z, "B", "N")
    
    return CRP2dcm(qcr)
end