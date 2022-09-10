function EP2dcm(q) where T<:Number
    return [q[4]^2 + q[1]^2 - q[2]^2 - q[3]^2   2*(q[1]*q[2] + q[4]*q[3])           2*(q[1]*q[3] - q[4]*q[2])
            2*(q[1]*q[2] - q[4]*q[3])           q[4]^2 - q[1]^2 + q[2]^2 - q[3]^2   2*(q[2]*q[3] + q[4]*q[1]) 
            2*(q[1]*q[3] + q[4]*q[2])           2*(q[2]*q[3] - q[4]*q[1])           q[4]^2 - q[1]^2 - q[2]^2 + q[3]^2]
end

function dcm2EP(dcm) #where T<:Number
    trace = dcm[1,1] + dcm[2,2] + dcm[3,3]
    q = sqrt.([ 1 - trace + 2*dcm[1,1]
                1 - trace + 2*dcm[2,2]
                1 - trace + 2*dcm[3,3]
                1 + trace]./4)
    index = [1 2 3 4][findfirst(q.==maximum(q))]
    if index == 4
        for j = 1:3
            i    = (1:3)[1:end .!= j]
            q[j] = (dcm[i[1],i[2]] - dcm[i[2],i[1]])/(4*q[index])*(-1).^(j + 1)
        end
    else
        for j = (1:3)[1:end .!= index]
            q[j] = (dcm[index,j] + dcm[j,index])/(4*q[index])
        end
        i    = (1:3)[1:end .!= index]
        q[4] = (dcm[i[1],i[2]] - dcm[i[2],i[1]])/(4*q[index])*(-1).^(index + 1)
    end
    return q
end

function dcm2EP_Mk2(dcm; η=0) #where T<:Number
    trace = dcm[1,1] + dcm[2,2] + dcm[3,3]
    check = [2*dcm[1,1] - trace
             2*dcm[2,2] - trace
             2*dcm[3,3] - trace
             trace]

    q = zeros(4,1)
    for i = 1:4
        k = (1:3)[1:end .!= i]
        sgn = sign((dcm[k[1],k[2]] - dcm[k[2],k[1]])*(-1).^(i + 1))*(i != 4) + (i == 4)
        if check[i] > η
            q[i] = 0.5*sqrt(1 + check[i])*sgn
        else
            q[i] = 0.5*sqrt(((dcm[3,2] + dcm[2,3]*(-1)^(i == 1 || i == 4))^2 + (dcm[1,2] + dcm[2,1]*(-1)^(i == 3 || i == 4))^2 + (dcm[3,1] + dcm[1,3]*(-1)^(i == 2 || i == 4))^2)/(3 - check[i]))*sgn
        end
    end
    return q
end

function EP2PRV()

end

function PRV2EP()

end

function SumEP(q1, q2)
    q = [q2[4]  q2[2] -q2[1] q2[3]
        -q2[2]  q2[4]  q2[3] q2[1]
         q2[1] -q2[3]  q2[4] q2[2]
        -q2[3] -q2[1] -q2[2] q2[4]]*q1
        return q
end

function EP2om(q) where T<:Number
    dq = inv([q[4] -q[2]  q[1] q[3]
              q[2]  q[4] -q[3] q[1]
             -q[1]  q[3]  q[4] q[2]
             -q[3] -q[1] -q[2] q[4]]./2)
    return dq
end

function om2EP(q) where T<:Number
    dq = [q[4] -q[2]  q[1]             # [q[4] -q[2]  q[1] q[3]
          q[2]  q[4] -q[3]             #  q[2]  q[4] -q[3] q[1]
         -q[1]  q[3]  q[4]             # -q[1]  q[3]  q[4] q[2]
         -q[3] -q[1] -q[2]]./2         # -q[3] -q[1] -q[2] q[4]]
    return dq
end

# if check[1] > η
#     q[1] = 0.5*sqrt(1 + check[1])*sign(dcm[3,2]-dcm[2,3])
# else
#     q[1] = 0.5*sqrt(((dcm[3,2] - dcm[2,3])^2 + (dcm[1,2] + dcm[2,1])^2 + (dcm[3,1] + dcm[1,3])^2)/(3 - check[1]))*sign(dcm[2,3] - dcm[3,2])
# end

# if check[2] > η
#     q[2] = 0.5*sqrt(1 + check[2])*sign(dcm[3,2]-dcm[2,3])
# else
#     q[2] = 0.5*sqrt(((dcm[3,2] + dcm[2,3])^2 + (dcm[1,2] + dcm[2,1])^2 + (dcm[3,1] - dcm[1,3])^2)/(3 - check[2]))*sign(dcm[3,1] - dcm[1,3])
# end

# if check[3] > η
#     q[3] = 0.5*sqrt(1 + check[3])*sign(dcm[3,2]-dcm[2,3])
# else
#     q[3] = 0.5*sqrt(((dcm[3,2] + dcm[2,3])^2 + (dcm[1,2] - dcm[2,1])^2 + (dcm[3,1] + dcm[1,3])^2)/(3 - check[3]))*sign(dcm[1,2] - dcm[2,1])
# end

# if check[4] > η
#     q[4] = 0.5*sqrt(1 + check[4])
# else
#     q[4] = 0.5*sqrt(((dcm[3,2] - dcm[2,3])^2 + (dcm[1,2] - dcm[2,1])^2 + (dcm[3,1] - dcm[1,3])^2)/(3 - check[4]))
# end
