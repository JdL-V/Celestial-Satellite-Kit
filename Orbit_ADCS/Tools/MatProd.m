function var = MatProd(A)
    N = size(A,2)
    Mat = eye(N,N)
    for k = 1:size(A,3)
        Mat = A(:,:,k)*Mat
    end
    var = Mat
end