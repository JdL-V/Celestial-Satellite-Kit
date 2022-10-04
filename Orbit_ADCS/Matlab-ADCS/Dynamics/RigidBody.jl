using LinearAlgebra
include("../Tools/MathTools.jl")
include("Continous.jl")

function getInertiaTensor(M::Vector{Float64}, R::Vector{Vector{Float64}})
    Io = zeros(3,3)
    for j::Int64 = 1:length(M)
        Io = Io - M[j].*cross2mat(R)^2
    end
    return I
end

function translateInertiaFrame(Ig::Matrix{Float64}, M::Vector{Float64}, Rp::Vector{Float64})
    return Ig + sum(M)*cross2mat(Rp)*cross2mat(Rp)'
end

function rotateInertiaFrame(Io::Matrix{Float64}, dcm::Matrix{Float64})
    return dcm*Io*dcm'
end

function getPrincipalInertia(Io::Matrix{Float64})
    dcm::Matrix{Float64} = eigvecs(Io)'
    Lam, dcm = eigenSort(eigvals(Io), dcm)
    # Right handed coordinate frame check
    dcm[3,:] = cross(dcm[1,:], dcm[2,:])
    return diagm(Lam), dcm
end

function getRotKineticE(Io::Matrix{Float64}, om::Vector{Float64})
    return 0.5*om'*Io*om
end