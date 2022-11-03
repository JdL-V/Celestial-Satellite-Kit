function var = getLinearMom(M, dR)
    dRc = getCM(M, dR);
    var = sum(M).*dRc;
end