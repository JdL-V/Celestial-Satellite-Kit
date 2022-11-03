function var = getKineticE(M, dR)
    dRc = getCM(M, dR);
    T1 = 0.5*sum(M)*dot(dRc, dRc);

    T2 = 0;
    T = 0;
    for j = 1:length(M)
        % T = T + 0.5*M(j)*dot(dR(j), dR(j))
        T2 = T2 + 0.5*M(j)*dot(dR(j) - dRc, dR(j) - dRc);
    end
    T = T1 + T2;
    var = [T, T1, T2];
end