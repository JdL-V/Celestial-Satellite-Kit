function EP2dcm(q::Vector{Float64})
    return [q[1]^2 + q[2]^2 - q[3]^2 - q[4]^2   2*(q[2]*q[3] + q[1]*q[4])           2*(q[2]*q[4] - q[1]*q[3])
            2*(q[2]*q[3] - q[1]*q[4])           q[1]^2 - q[2]^2 + q[3]^2 - q[4]^2   2*(q[3]*q[4] + q[1]*q[2]) 
            2*(q[2]*q[4] + q[1]*q[3])           2*(q[3]*q[4] - q[1]*q[2])           q[1]^2 - q[2]^2 - q[3]^2 + q[4]^2]
end

function sheppard(dcm::Matrix{Float64})
    trace::Float64 = dcm[1,1] + dcm[2,2] + dcm[3,3]
    q::Vector{Float64}= sqrt.([ 1 + trace
                                1 - trace + 2*dcm[1,1]
                                1 - trace + 2*dcm[2,2]
                                1 - trace + 2*dcm[3,3]]./4)
    index::Int8 = [1 2 3 4][findfirst(q.==maximum(q))]
    if index == 1
        for j::Int8 = 1:3
            i::Vector{Int8} = (1:3)[1:end .!= j]
            q[j+1]  = (dcm[i[1],i[2]] - dcm[i[2],i[1]])/(4*q[index])*(-1).^(j + 1)
        end
    else
        i = (1:3)[1:end .!= index-1]
        q[1] = (dcm[i[1],i[2]] - dcm[i[2],i[1]])/(4*q[index])*(-1).^index
        for j::Int8 = (1:3)[1:end .!= index-1]
            q[j+1] = (dcm[index-1,j] + dcm[j,index-1])/(4*q[index])
        end
        q = q.*sign(q[1])
    end
    return q
end

function dcm2EP(dcm::Matrix{Float64}; η::Float64=0.)
    trace::Float64 = dcm[1,1] + dcm[2,2] + dcm[3,3]
    check::Vector{Float64} = [trace,
                              2*dcm[1,1] - trace,
                              2*dcm[2,2] - trace,
                              2*dcm[3,3] - trace]

    q::Vector{Float64} = zeros(4)
    for i::Int8 = 1:4
        k::Vector{Int8} = (1:3)[1:end .!= (i-1)]
        sgn::Int8 = sign((dcm[k[1],k[2]] - dcm[k[2],k[1]])*(-1).^i)*(i != 1) + (i == 1)

        if check[i] > η
            q[i] = 0.5*sqrt(1 + check[i])*sgn
        else
            q[i] = 0.5*sqrt((   (dcm[3,2] + dcm[2,3]*(-1)^(i == 2 || i == 1))^2 
                              + (dcm[1,2] + dcm[2,1]*(-1)^(i == 4 || i == 1))^2 
                              + (dcm[3,1] + dcm[1,3]*(-1)^(i == 3 || i == 1))^2) / (3 - check[i]))*sgn
        end
    end
    return q
end

function EP2PRV(q::Vector{Float64})
    ϕ::Float64 = 2*acos(q[1])
    e::Float64 = q[2:4]./sin(ϕ/2)
    return ϕ, e
end

function PRV2EP(ϕ::Float64, e::Vector{Float64})
    return [cos(ϕ/2), e[1]sin(ϕ/2), e[2]sin(ϕ/2), e[3]sin(ϕ/2)]
end

function SumEP(q1::Vector{Float64}, q2::Vector{Float64})
    return [q2[1] -q2[2] -q2[3] -q2[4]
            q2[2]  q2[1]  q2[4] -q2[3]
            q2[3] -q2[4]  q2[1]  q2[2]
            q2[4]  q2[3] -q2[2]  q2[1]]*q1
end

function EP2om(q::Vector{Float64})
    return inv([q[1] -q[2] -q[3] -q[4]
                q[2]  q[1] -q[4]  q[3]
                q[3]  q[4]  q[1] -q[2]
                q[4] -q[3]  q[2]  q[1]]./2)
end

function om2EP(q::Vector{Float64})
    return [-q[2] -q[3] -q[4]             # [q[1] -q[2] -q[3] -q[4]
             q[1] -q[4]  q[3]             #  q[2]  q[1] -q[4]  q[3]
             q[4]  q[1] -q[2]             #  q[3]  q[4]  q[1] -q[2]
            -q[3]  q[2]  q[1]]./2         #  q[4] -q[3]  q[2]  q[1]]
end

function EPdiff(nn::Int64, tmax::Float64, u0::Vector{Float64}, f::Function, sch::Function)
    h::Float64  = tmax/(nn-1)
    th::LinRange{Float64, Int64} = LinRange(0, tmax, nn)
    u::Matrix{Float64}  = zeros(nn, length(u0))
    u[1,:] = u0

    for n::Int64 = 1:(nn-1)
        u[n+1,:] = normalize(sch(h, u[n,:], th[n], f))
    end
    return th, u
end