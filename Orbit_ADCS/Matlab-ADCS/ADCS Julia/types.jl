mutable struct DCM
    Mat::Matrix{Float64}
    outFrame::String
    inFrame::String
end

mutable struct EulerAng
    Angles::Vector{Float64}
    Sequence::Vector{Int64}
    outFrame::String
    inFrame::String
end

mutable struct PRV
    Angle::Float64
    x::Vector{Float64}
    outFrame::String
    inFrame::String
end

mutable struct quaternion
    x::Vector{Float64}
    outFrame::String
    inFrame::String
end

mutable struct CRP
    x::Vector{Float64}
    outFrame::String
    inFrame::String
end

mutable struct MRP
    x::Vector{Float64}
    outFrame::String
    inFrame::String
end

mutable struct FMatrix
    Mat::Matrix{Float64}
    Frame::String
end

mutable struct FVector
    Vec::Vector{Float64}
    Frame::String
end

include("Tools/MathTools.jl")

function dcmmul(dcm2::DCM, dcm1::DCM)
    if dcm1.outFrame != dcm2.inFrame
        error("dcms are not frame-compatible")
    end
    dcm = DCM(dcm2.Mat*dcm1.Mat, dcm2.outFrame, dcm1.inFrame)
    return dcm
end

function dcmtrs(dcm::DCM)
    dcmout = DCM(trs(dcm.Mat), dcm.inFrame, dcm.outFrame)
    return dcmout
end