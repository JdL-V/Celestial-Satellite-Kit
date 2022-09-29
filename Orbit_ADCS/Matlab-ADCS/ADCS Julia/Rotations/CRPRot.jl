include("../Tools/MathTools.jl")
include("EPRot.jl")

function CRP2dcm(qcr::CRP)
    dcm = DCM([ 1 + qcr.x[1]^2 - qcr.x[2]^2 - qcr.x[3]^2    2*(qcr.x[1]*qcr.x[2] + qcr.x[3])            2*(qcr.x[1]*qcr.x[3] - qcr.x[2])
                2*(qcr.x[1]*qcr.x[2] - qcr.x[3])            1 - qcr.x[1]^2 + qcr.x[2]^2 - qcr.x[3]^2    2*(qcr.x[2]*qcr.x[3] + qcr.x[1]) 
                2*(qcr.x[1]*qcr.x[3] + qcr.x[2])            2*(qcr.x[2]*qcr.x[3] - qcr.x[1])            1 - qcr.x[1]^2 - qcr.x[2]^2 + qcr.x[3]^2]./(1 + dot(qcr.x,qcr.x)),
            qcr.outFrame,
            qcr.inFrame)
    return dcm
    # return ((1 - dot(qcr,qcr))*Matrix(I, 3, 3) + 2*qcr*qcr' - 2*cross2mat(qcr))./(1 + dot(qcr,qcr))
end

function dcm2CRP(dcm::DCM)
    # ζ = sqrt(1 + dcm[1,1] + dcm[2,2] + dcm[3,3])
    # return [dcm[2,3] - dcm[3,2]
    #         dcm[3,1] - dcm[1,3] 
    #         dcm[1,2] - dcm[2,1]]./ζ^2
    return EP2CRP(dcm2EP(dcm))
end

function CRP2EP(qcr::CRP)
    q = quaternion([1, qcr.x[1:3]]./sqrt(1 + dot(qcr,qcr)),
                    qcr.outFrame,
                    qcr.inFrame)
    return q
end

function EP2CRP(q::quaternion)
    qcr = CRP(q.x[2:4]./q.x[1], 
            q.outFrame, 
            q.inFrame)
    return qcr
end

function CRP2PRV(qcr::CRP)
    T::Float64 = atan.(qcr.x)
    e::Vector{Float64} = normalize(T)
    prv = PRV(2*T[1]/e[1], e, qcr.outFrame, qcr.inFrame)
    return prv
end

function PRV2CRP(prv::PRV)
    qcr = CRP(tan(prv.Angle/2).*prv.x, prv.outFrame, prv.inFrame)
    return qcr
end

function SumCRP(qcr1::CRP, qcr2::CRP)
    qcr = CRP((qcr2 + qcr1 - cross(qcr2,qcr1))./(1 - dot(qcr2,qcr1)),
                qcr2.outFrame,
                qcr1.inFrame)
    return qcr
end

function CRP2om(qcr::Vector{Float64})
    return 2/(1 + dot(qcr, qcr)).*(Matrix(I, 3, 3) - cross2mat(qcr))
end

function om2CRP(qcr::Vector{Float64})
    return [1 + qcr[1]^2            qcr[1]*qcr[2] - qcr[3]  qcr[1]*qcr[3] + qcr[2]
            qcr[2]*qcr[1] + qcr[3]  1 + qcr[2]^2            qcr[2]*qcr[3] - qcr[1]
            qcr[3]*qcr[1] - qcr[2]  qcr[3]*qcr[2] + qcr[1]  1 + qcr[3]^2          ].*0.5
    # return 0.5.*(Matrix(I, 3, 3) + cross2mat(qcr) + qcr*qcr')
end