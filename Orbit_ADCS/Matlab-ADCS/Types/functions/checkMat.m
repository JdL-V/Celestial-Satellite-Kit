function var = checkMat(Mat, Dim)
    if ~(prod(size(Mat) == Dim))
        error('Matrix dimesions not %d by %d', Dim(1), Dim(2))
    end
end