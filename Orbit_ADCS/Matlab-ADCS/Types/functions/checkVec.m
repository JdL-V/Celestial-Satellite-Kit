function var = checkVec(Vec, Dim)
    if ~(length(Vec(:)) == Dim)
        error('Vector dimesion not %d', Dim)
    end
end