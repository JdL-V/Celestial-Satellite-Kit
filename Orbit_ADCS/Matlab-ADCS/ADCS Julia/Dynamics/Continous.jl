using LinearAlgebra
function getCM(M::Vector{Float64}, R::Vector{Vector{Float64}})
    # Works for position, speed and acceleration
    CM::Vector{Float64} = zeros(3)
    for j = 1:length(M)
        CM = CM + M[j].*R[j]
    end
    return CM/sum(M)
end

function getForce(M::Vector{Float64}, ddRc::Vector{Vector{Float64}})
    return sum(M).*ddRc
end

function getKineticE(M::Vector{Float64}, dR::Vector{Vector{Float64}})
    dRc::Vector{Float64} = getCM(M, dR)
    T1::Float64 = 0.5*sum(M)*dot(dRc, dRc)

    T2::Float64 = 0
    T::Float64 = 0
    for j::Int64 = 1:length(M)
        # T = T + 0.5*M[j]*dot(dR[j], dR[j])
        T2 = T2 + 0.5*M[j]*dot(dR[j] - dRc, dR[j] - dRc)
    end
    T = T1 + T2
    return T, T1, T2
end

function getLinearMom(M::Vector{Float64}, dR::Vector{Vector{Float64}})
    dRc::Vector{Float64} = getCM(M, dR)
    return sum(M).*dRc
end

function getAngularMom(M::Vector{Float64}, R::Vector{Vector{Float64}}, dR::Vector{Vector{Float64}}, Rp::Vector{Float64}, dRp::Vector{Float64})
    H::Vector{Float64} = zeros(3)
    for j::Int64 = 1:length(M)
        H = H + M[j].*cross(R[j] - Rp, dR[j] - dRp)
    end
    return H
end

function getTorque(M::Vector{Float64}, R::Vector{Vector{Float64}}, Rp::Vector{Float64}, F::Vector{Float64})
    Lp::Vector{Float64} = zeros(3)
    for j::Int64 = 1:length(M)
        Lp = Lp + cross(R[j] - Rp, F)
    end
    return Lp
end

# function getInertialtDev(M::Vector{Float64}, R::Vector{Vector{Float64}}, Rp::Vector{Float64}, F::Vector{Float64})
#     Lp::Vector{Float64} = zeros(3)
#     for j::Int64 = 1:length(M)
#         Lp = Lp + cross(R[j] - Rp, F)
#     end
#     return Lp
# end