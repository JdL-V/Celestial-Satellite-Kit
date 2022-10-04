function var = qMethod(vkb, vkn, weight)
    TP = types;
    EP = EPRot;
    B = zeros(3,3);
    for i = 1:size(vkb, 2)
        B = B + weight(i)*vkb(:,i)*vkn(:,i)'
    end

    S = B + B';
    sig = trace(B);
    Z = [B(2,3) - B(3,2)
         B(3,1) - B(1,3)
         B(1,2) - B(2,1)]

    K = [sig Z'
         Z   S-sig.*eye(3,3)];
    
    lambda = eigvals(K);
    max_lambda = findfirst(real(lambda) == maximum(real(lambda)));
    q = TP.quaternion(real(eigvecs(K)(:,max_lambda)), "B", "N");
    
    var = EP.EP2dcm(q);
end