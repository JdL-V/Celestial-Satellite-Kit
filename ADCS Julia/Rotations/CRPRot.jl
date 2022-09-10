include("../Tools/MathTools.jl")
include("EPRot.jl")

function CRP2dcm(qcr::Vector{Float64})
    return [1 + qcr[1]^2 - qcr[2]^2 - qcr[3]^2    2*(qcr[1]*qcr[2] + qcr[3])            2*(qcr[1]*qcr[3] - qcr[2])
            2*(qcr[1]*qcr[2] - qcr[3])            1 - qcr[1]^2 + qcr[2]^2 - qcr[3]^2    2*(qcr[2]*qcr[3] + qcr[1]) 
            2*(qcr[1]*qcr[3] + qcr[2])            2*(qcr[2]*qcr[3] - qcr[1])            1 - qcr[1]^2 - qcr[2]^2 + qcr[3]^2]./(1 + dot(qcr,qcr))
    # return ((1 - dot(qcr,qcr))*Matrix(I, 3, 3) + 2*qcr*qcr' - 2*cross2mat(qcr))./(1 + dot(qcr,qcr))
end

function dcm2CRP(dcm::Matrix{Float64})
    # ζ = sqrt(1 + dcm[1,1] + dcm[2,2] + dcm[3,3])
    # return [dcm[2,3] - dcm[3,2]
    #         dcm[3,1] - dcm[1,3] 
    #         dcm[1,2] - dcm[2,1]]./ζ^2
    return EP2CRP(dcm2EP(dcm))
end

function CRP2EP(qcr::Vector{Float64})
    return [1, qcr[1:3]]./sqrt(1 + dot(qcr,qcr))
end

function EP2CRP(q::Vector{Float64})
    return q[2:4]./q[1]
end

function CRP2PRV(qcr::Vector{Float64})
    T::Float64 = atan.(qcr)
    e::Vector{Float64} = normalize(T)
    return 2*T[1]/e[1], e
end

function PRV2CRP(ϕ::Float64, e::Vector{Float64})
    return tan(ϕ/2).*e
end

function SumCRP(qcr1::Vector{Float64}, qcr2::Vector{Float64})
    return (qcr2 + qcr1 - cross(qcr2,qcr1))./(1 - dot(qcr2,qcr1))
end

function CRP2om(qcr::Vector{Float64})
    return 2/(1 + dot(qcr, qcr)).*(Matrix(I, 3, 3) - cross2mat(qcr))
end

function om2CRP(qcr::Vector{Float64})
    return [1 + qcr[1]^2                qcr[1]*qcr[2] - qcr[3]      qcr[1]*qcr[3] + qcr[2]
            qcr[2]*qcr[1] + qcr[3]      1 + qcr[2]^2                qcr[2]*qcr[3] - qcr[1]
            qcr[3]*qcr[1] - qcr[2]      qcr[3]*qcr[2] + qcr[1]      1 + qcr[3]^2            ].*0.5
    # return 0.5.*(Matrix(I, 3, 3) + cross2mat(qcr) + qcr*qcr')
end