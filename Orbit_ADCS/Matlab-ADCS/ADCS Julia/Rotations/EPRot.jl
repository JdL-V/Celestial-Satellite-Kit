include("../types.jl")

function EP2dcm(q::quaternion)
    dcm = DCM([q.x[1]^2 + q.x[2]^2 - q.x[3]^2 - q.x[4]^2   2*(q.x[2]*q.x[3] + q.x[1]*q.x[4])           2*(q.x[2]*q.x[4] - q.x[1]*q.x[3])
            2*(q.x[2]*q.x[3] - q.x[1]*q.x[4])           q.x[1]^2 - q.x[2]^2 + q.x[3]^2 - q.x[4]^2   2*(q.x[3]*q.x[4] + q.x[1]*q.x[2]) 
            2*(q.x[2]*q.x[4] + q.x[1]*q.x[3])           2*(q.x[3]*q.x[4] - q.x[1]*q.x[2])           q.x[1]^2 - q.x[2]^2 - q.x[3]^2 + q.x[4]^2],
            q.outFrame,
            q.inFrame)
    return dcm
end

function sheppard(dcm::DCM)
    trace::Float64 = dcm.Mat[1,1] + dcm.Mat[2,2] + dcm.Mat[3,3]
    q = quaternion( sqrt.([ 1 + trace
                            1 - trace + 2*dcm.Mat[1,1]
                            1 - trace + 2*dcm.Mat[2,2]
                            1 - trace + 2*dcm.Mat[3,3]]./4),
                    dcm.outFrame,
                    dcm.inFrame)

    index::Int8 = [1 2 3 4][findfirst(q.x .== maximum(q.x))]
    if index == 1
        for j::Int8 = 1:3
            i::Vector{Int8} = (1:3)[1:end .!= j]
            q.x[j+1]  = (dcm.Mat[i[1],i[2]] - dcm.Mat[i[2],i[1]])/(4*q.x[index])*(-1).^(j + 1)
        end
    else
        i = (1:3)[1:end .!= index-1]
        q.x[1] = (dcm.Mat[i[1],i[2]] - dcm.Mat[i[2],i[1]])/(4*q.x[index])*(-1).^index
        for j::Int8 = (1:3)[1:end .!= index-1]
            q.x[j+1] = (dcm.Mat[index-1,j] + dcm.Mat[j,index-1])/(4*q.x[index])
        end
        q.x = (q.x).*sign(q.x[1])
    end
    return q
end

function dcm2EP(dcm::DCM; η::Float64=0.)
    trace::Float64 = dcm.Mat[1,1] + dcm.Mat[2,2] + dcm.Mat[3,3]
    check::Vector{Float64} = [trace,
                              2*dcm.Mat[1,1] - trace,
                              2*dcm.Mat[2,2] - trace,
                              2*dcm.Mat[3,3] - trace]

    q = quaternion(zeros(4), dcm.outFrame, dcm.inFrame)
    for i::Int8 = 1:4
        k::Vector{Int8} = (1:3)[1:end .!= (i-1)]
        sgn::Int8 = sign((dcm.Mat[k[1],k[2]] - dcm.Mat[k[2],k[1]])*(-1).^i)*(i != 1) + (i == 1)

        if check[i] > η
            q.x[i] = 0.5*sqrt(1 + check[i])*sgn
        else
            q.x[i] = 0.5*sqrt(( (dcm.Mat[3,2] + dcm.Mat[2,3]*(-1)^(i == 2 || i == 1))^2 
                              + (dcm.Mat[1,2] + dcm.Mat[2,1]*(-1)^(i == 4 || i == 1))^2 
                              + (dcm.Mat[3,1] + dcm.Mat[1,3]*(-1)^(i == 3 || i == 1))^2) / (3 - check[i]))*sgn
        end
    end
    return q
end

function EP2PRV(q::quaternion)
    angle::Float64 = 2*acos(q.x[1])
    prv = PRV(
        angle,
        q.x[2:4]./sin(angle/2),
        q.outFrame,
        q.inFrame)
    return prv
end

function PRV2EP(prv::PRV)
    q = quaternion(
        [cos(prv.Angle/2), prv.x[1]sin(prv.Angle/2), prv.x[2]sin(prv.Angle/2), prv.x[3]sin(prv.Angle/2)],
        prv.outFrame,
        prv.inFrame)
    return q
end

function SumEP(q1::quaternion, q2::quaternion)
    q = quaternion([q2[1] -q2[2] -q2[3] -q2[4]
                    q2[2]  q2[1]  q2[4] -q2[3]
                    q2[3] -q2[4]  q2[1]  q2[2]
                    q2[4]  q2[3] -q2[2]  q2[1]]*q1,
                    q2.outFrame,
                    q1.inFrame)
    return q
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
        u[n+1,:] = sch(h, u[n,:], th[n], f)
        # (1:3 to use the equations with omega rates)
        u[n+1,1:4] = normalize(u[n+1,1:4]) 
    end
    return th, u
end