function var = om2CRP(qcr)
    var = [1 + qcr(1)^2            qcr(1)*qcr(2) - qcr(3)  qcr(1)*qcr(3) + qcr(2)
           qcr(2)*qcr(1) + qcr(3)  1 + qcr(2)^2            qcr(2)*qcr(3) - qcr(1)
           qcr(3)*qcr(1) - qcr(2)  qcr(3)*qcr(2) + qcr(1)  1 + qcr(3)^2          ].*0.5;
    % var = 0.5.*(Matrix(I, 3, 3) + skewsym(qcr) + qcr*qcr')
end