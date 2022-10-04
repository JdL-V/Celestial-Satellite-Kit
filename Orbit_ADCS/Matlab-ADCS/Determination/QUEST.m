function var = QUEST(vkb, vkn, weight)
    TP = types;
    CRP = CRPRot;
    B = zeros(3,3)
    for i = 1:size(vkb, 2)
        B = B + weight(i)*vkb(:,i)*vkn(:,i)';
    end

    S = B + B';
    sig = trace(B);
    Z = [B(2,3) - B(3,2)
         B(3,1) - B(1,3)
         B(1,2) - B(2,1)];

    K = [sig Z'
         Z   S-sig.*eye(3,3)];
    
    f(s) = det(K - s(1).*eye(4,4));
    lambda0 = sum(weight);
    lambdaopt = newton(lambda0, f);

    qcr = TP.CRP(inv((lambda_opt + sig)*eye(3,3) - S)*Z, "B", "N");
    
    var = CRP.CRP2dcm(qcr);
end