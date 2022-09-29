include("../Tools/MathTools.jl")
include("../types.jl")

function PRV2dcm(prv::PRV)
    Σ::Float64 = 1 - cos(prv.Angle)
    dcm = DCM(zeros(3,3), prv.outFrame, prv.inFrame)
    # needs revision
    for i::Int8 = 1:3
        for j::Int8 = 1:3
            dcm.Mat[i,j] = prv.x[i]*prv.x[j]*Σ
            if i == j
                dcm.Mat[i,j] = dcm.Mat[i,j] + cos(prv.Angle)
            else
                dcm.Mat[i,j] = dcm.Mat[i,j] + (-1)^(i + j + (i < j))*sin(prv.Angle)*prv.x[([1 2 3]*prod.(eachcol([1 2 3].!=[i; j])))[1]]
            end
        end
    end
    return dcm
end

function dcm2PRV(dcm::DCM)
    prv = PRV(0., [0.,0.,0.], dcm.outFrame, dcm.inFrame)
    prv.Angle::Float64 = acos(0.5*(dcm.Mat[1,1] + dcm.Mat[2,2] + dcm.Mat[3,3] - 1))
    prv.x::Vector{Float64} = 1/(2*sin(prv.Angle))*[dcm.Mat[2,3]-dcm.Mat[3,2], dcm.Mat[3,1]-dcm.Mat[1,3], dcm.Mat[1,2]-dcm.Mat[2,1]]
    return prv
end

function PRV2om(prv::PRV)
    τ = cross2mat(prv.x*prv.Angle)
    return (Matrix(I, 3, 3) + 0.5.*τ + 1/prv.Angle^2*(1 - prv.Angle/2*cot(prv.Angle/2)).*τ^2)
end

function om2PRV(prv::PRV)
    τ::Matrix{Float64} = cross2mat(prv.x.*prv.Angle)
    return (Matrix(I, 3, 3) - ((1 - cos(prv.Angle))/prv.Angle^2).*τ + ((prv.Angle - sin(prv.Angle))/prv.Angle^3).*τ^2)
end