function var = OLAE2(vkb, vkn, weight)
    N = size(vkb, 2);
    weight = weight/sum(weight(1:N));

    B = zeros(3,3);
    Mm = zeros(3,3);   
    for i = 1:N
        B = B + weight(i)*vkb(:,i)*vkn(:,i)';
        Mm = Mm + weight(i)*(vkb(:,i)*vkb(:,i)' + vkn(:,i)*vkn(:,i)');
    end

    S = B + B';
    sig = trace(B);
    Z = [B(2,3) - B(3,2)
         B(3,1) - B(1,3)
         B(1,2) - B(2,1)]; 

    Mw = 0.5.*S - (sig + 1).*eye(3,3);
    M = (0.5.*Mm + Mw)

    g = M\-Z
    
    qcr = CRP(g, "B", "N");
    var = CRP2dcm(qcr);
end