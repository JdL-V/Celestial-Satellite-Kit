function var = getAngularMom(M, R, dR, Rp, dRp)
    H = zeros(3);
    for j = 1:length(M)
        H = H + M(j).*cross(R(j) - Rp, dR(j) - dRp);
    end
    var = H;
end