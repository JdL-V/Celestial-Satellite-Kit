include("../Tools/MathTools.jl")

function PRV2dcm(ϕ::Float64, e::Vector{Float64})
    Σ::Float64 = 1 - cos(ϕ)
    dcm::Matrix{Float64} = zeros(3,3)
    # needs revision
    for i::Int8 = 1:3
        for j::Int8 = 1:3
            dcm[i,j] = e[i]*e[j]*Σ
            if i == j
                dcm[i,j] = dcm[i,j] + cos(ϕ)
            else
                dcm[i,j] = dcm[i,j] + (-1)^(i + j + (i < j))*sin(ϕ)*e[([1 2 3]*prod.(eachcol([1 2 3].!=[i; j])))[1]]
            end
        end
    end
    return dcm
end

function dcm2PRV(dcm::Matrix{Float64})
    ϕ::Float64 = acos(0.5*(dcm[1,1] + dcm[2,2] + dcm[3,3] - 1))
    e::Vector{Float64} = 1/(2*sin(ϕ))*[dcm[2,3]-dcm[3,2], dcm[3,1]-dcm[1,3], dcm[1,2]-dcm[2,1]]
    return ϕ, e
end

function PRV2om(ϕ::Float64, e::Vector{Float64})
    τ = cross2mat(e*ϕ)
    return (Matrix(I, 3, 3) + 0.5.*τ + 1/ϕ^2*(1 - ϕ/2*cot(ϕ/2)).*τ^2)
end

function om2PRV(ϕ::Float64, e::Vector{Float64})
    τ::Matrix{Float64} = cross2mat(e.*ϕ)
    return (Matrix(I, 3, 3) - ((1 - cos(ϕ))/ϕ^2).*τ + ((ϕ - sin(ϕ))/ϕ^3).*τ^2)
end