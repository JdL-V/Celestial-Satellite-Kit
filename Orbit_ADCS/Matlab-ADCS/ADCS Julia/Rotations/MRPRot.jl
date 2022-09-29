include("../Tools/MathTools.jl")
include("EPRot.jl")
include("../types.jl")

function MRP2dcm(qmr::MRP)
    qmrMat::Matrix{Float64} = cross2mat(qmr.x)
    qmrSq::Float64 = dot(qmr.x,qmr.x)
    dcm = DCM(Matrix(I, 3, 3) + (8*qmrMat^2 - 4*(1 - qmrSq)*qmrMat)/(1 + qmrSq)^2,
            qmr.outFrame,
            qmr.inFrame)
    return dcm
end

function dcm2MRP(dcm::DCM)
    # ζ = sqrt(1 + dcm[1,1] + dcm[2,2] + dcm[3,3])
    # return [dcm[2,3] - dcm[3,2]
    #         dcm[3,1] - dcm[1,3] 
    #         dcm[1,2] - dcm[2,1]]./(ζ*(ζ + 2))
    return EP2MRP(dcm2EP(dcm))
end

function MRP2EP(qmr::MRP)
    Σ::Float64 = dot(qmr.x,qmr.x)
    q = quaternion([(1 - Σ)/(1 + Σ) 2*qmr[1:3]'/(1 + Σ)],
                    qmr.outFrame,
                    qmr.inFrame)
    return q
end

function EP2MRP(q::quaternion)
    qmr = MRP(q.x[2:4]./(1 + q.x[1]),
                q.outFrame,
                q.inFrame)
    return qmr
end

function MRP2CRP(qmr::MRP)
    qcr = CRP(qmr.x./(1 - dot(qmr.x,qmr.x)),
                qmr.outFrame,
                qmr.inFrame)
    return qcr
end

function CRP2MRP(qcr::CRP)
    qmr = MRP(qcr./(1 + sqrt(1 + dot(qcr,qcr))),
                qcr.outFrame,
                qcr.inFrame)
    return qmr
end

function MRP2PRV(qmr::MRP)
    T::Vector{Float64} = atan.(qmr)
    e::Vector{Float64} = normalize(T)
    prv = PRV(4*T[1]/e[1], e,
                qmr.outFrame,
                qmr.inFrame)
    return prv
end

function PRV2MRP(prv::PRV)
    qmr = MRP(tan(prv.Angle/4).*prv.x,
                prv.outFrame,
                prv.inFrame)
    return qmr
end

function SumMRP(qmr1::MRP, qmr2::MRP)
    qmr1Sq = dot(qmr1.x,qmr1.x)
    qmr2Sq = dot(qmr2.x,qmr2.x)
    qmr = MRP(((1 - qmr1Sq)*qmr2.x + (1 - qmr2Sq)*qmr1.x - 2*cross(qmr2.x,qmr1.x))
            ./(1 + qmr1Sq*qmr2Sq - 2*dot(qmr1.x,qmr2.x)),
            qmr2.outFrame,
            qmr1.inFrame)
    return qmr
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
        # Avoid singularity (1:3 in order to use the equations with omega rates)
        if dot(u[n+1,1:3], u[n+1,1:3]) > 1
            u[n+1,1:3] = -u[n+1,1:3]./dot(u[n+1,1:3],u[n+1,1:3])
        end
    end
    return th, u
end