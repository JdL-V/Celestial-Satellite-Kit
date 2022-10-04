function var = OLAE(vkb, vkn, weight)
    TP = types;
    CRP = CRPRot;
    d = vec(vkb - vkn);
    N = length(d);
    S = zeros(N,3);
    for i = 1:Int(N/3)
        S(i:i+2,:) = cross2mat(vkb(:,i) + vkn(:,i));
    end
    W = eye(N,N)
    qcr = TP.CRP(inv(S'*W*S)*S'*W*d, "B", "N");

    var = CRP.CRP2dcm(qcr);
end