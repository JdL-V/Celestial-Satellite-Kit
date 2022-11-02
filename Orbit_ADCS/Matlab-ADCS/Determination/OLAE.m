function var = OLAE(vkb, vkn, weight)
    TP = types;
    CRP = CRPRot;
    d0 = vkb - vkn;
    d = reshape(d0,[numel(d0),1]);
    N = length(d);
    S = zeros(N,3);
    for i = 1:(N/3)
        S(i:i+2,:) = skewsym(vkb(:,i) + vkn(:,i));
        for j = 3*i - [2 1 0]
            W(j,j) = weight(i);
        end
    end
    SWS = (S'*W*S);
    if 1/cond(SWS) < 1e-4
        warning("Matrix is close to singular or badly scaled. Results may be inaccurate. Changing to QUEST method")
        var = QUEST(vkb, vkn, weight);
        return
    end
    qcr = TP.CRP(SWS\(S'*W*d), "B", "N");

    var = CRP.CRP2dcm(qcr);
end