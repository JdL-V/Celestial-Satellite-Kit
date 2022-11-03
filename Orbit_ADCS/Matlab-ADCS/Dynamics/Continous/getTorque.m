function var = getTorque(M, R, Rp, F)
    Lp = zeros(3,1);
    for j = 1:length(M)
        Lp = Lp + cross(R(j) - Rp, F);
    end
    var = Lp;
end