function trs(A::Matrix)
    return Matrix(A')
end

function eigenSort(Vals::Vector{Float64}, Vecs::Matrix{Float64})
    Vecs = Vecs[sortperm(Vals, rev=true),:]
    Vals = sort(Vals, rev=true)
    return Vals, Vecs
end

function MatProd(A)
    N::Int64 = size(A,2)
    Mat::Matrix{Number} = Matrix(I,N,N)
    for k = 1:size(A,3)
        Mat = A[:,:,k]*Mat
    end
    return Mat
end