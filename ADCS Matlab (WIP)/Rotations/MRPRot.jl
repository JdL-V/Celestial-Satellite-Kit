include("../Tools/MathTools.jl")
include("EPRot.jl")

function MRP2dcm(qmr::Vector{Float64})
    qmrMat::Matrix{Float64} = cross2mat(qmr)
    qmrSq::Float64 = dot(qmr,qmr)
    return Matrix(I, 3, 3) + (8*qmrMat^2 - 4*(1 - qmrSq)*qmrMat)/(1 + qmrSq)^2
end

function dcm2MRP(dcm::Matrix{Float64})
    # ζ = sqrt(1 + dcm[1,1] + dcm[2,2] + dcm[3,3])
    # return [dcm[2,3] - dcm[3,2]
    #         dcm[3,1] - dcm[1,3] 
    #         dcm[1,2] - dcm[2,1]]./(ζ*(ζ + 2))
    return EP2MRP(dcm2EP(dcm))
end

function MRP2EP(qmr::Vector{Float64})
    Σ::Float64 = dot(qmr,qmr)
    return [(1 - Σ)/(1 + Σ) 2*qmr[1:3]'/(1 + Σ)]
end

function EP2MRP(q::Vector{Float64})
    return q[2:4]./(1 + q[1])
end

function MRP2CRP(qmr::Vector{Float64})
    return qmr./(1 - dot(qmr,qmr))
end

function CRP2MRP(qcr::Vector{Float64})
    return qcr./(1 + sqrt(1 + dot(qcr,qcr)))
end

function MRP2PRV(qmr::Vector{Float64})
    T::Vector{Float64} = atan.(qmr)
    e::Vector{Float64} = normalize(T)
    return 4*T[1]/e[1], e
end

function PRV2MRP(ϕ::Float64, e::Vector{Float64})
    return tan(ϕ/4).*e
end

function SumMRP(qmr1::Vector{Float64}, qmr2::Vector{Float64})
    qmr1Sq = dot(qmr1,qmr1)
    qmr2Sq = dot(qmr2,qmr2)
    return (((1 - qmr1Sq)*qmr2 + (1 - qmr2Sq)*qmr1 - 2*cross(qmr2,qmr1))
            ./(1 + qmr1Sq*qmr2Sq - 2*dot(qmr1,qmr2)))
end

function MRP2om(qmr::Vector{Float64})
    return (4/(1 + dot(qmr,qmr))^2)*(om2MRP(qmr)')
end

function om2MRP(qmr::Vector{Float64})
    qmrSq::Float64 = dot(qmr,qmr)
    return [1 - qmrSq + 2*qmr[1]^2      2*(qmr[1]*qmr[2] - qmr[3])  2*(qmr[1]*qmr[3] + qmr[2])
            2*(qmr[2]*qmr[1] + qmr[3])  1 - qmrSq + 2*qmr[2]^2      2*(qmr[2]*qmr[3] - qmr[1])
            2*(qmr[3]*qmr[1] - qmr[2])  2*(qmr[3]*qmr[2] + qmr[1])  1 - qmrSq + 2*qmr[3]^2    ]./4
    # return ((1 - qmrSq).*Matrix(I,3,3) + 2*cross2mat(qmr) + 2*qmr*qmr')./4
end

function MRPdiff(nn::Int64, tmax::Float64, u0::Vector{Float64}, f::Function, sch::Function)
    h::Float64  = tmax/(nn-1)
    th::LinRange{Float64, Int64} = LinRange(0, tmax, nn)
    u::Matrix{Float64}  = zeros(nn, length(u0))
    u[1,:] = u0

    for n::Int64 = 1:(nn-1)
        u[n+1,:] = sch(h, u[n,:], th[n], f)
        # Avoid singularity
        if dot(u[n+1,:], u[n+1,:]) > 1
            u[n+1,:] = -u[n+1,:]./dot(u[n+1,:],u[n+1,:])
        end
    end
    return th, u
end