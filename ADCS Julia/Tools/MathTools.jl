function cross2mat(w::Vector{Float64})
    return [0      -w[3]    w[2]
            w[3]    0      -w[1]
           -w[2]    w[1]    0]
end