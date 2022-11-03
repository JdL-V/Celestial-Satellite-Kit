function var = om2MRP(qmr)
    qmrSq = dot(qmr,qmr);
    var = [1 - qmrSq + 2*qmr(1)^2      2*(qmr(1)*qmr(2) - qmr(3))  2*(qmr(1)*qmr(3) + qmr(2))
           2*(qmr(2)*qmr(1) + qmr(3))  1 - qmrSq + 2*qmr(2)^2      2*(qmr(2)*qmr(3) - qmr(1))
           2*(qmr(3)*qmr(1) - qmr(2))  2*(qmr(3)*qmr(2) + qmr(1))  1 - qmrSq + 2*qmr(3)^2    ]./4;
    % var = ((1 - qmrSq).*Matrix(I,3,3) + 2*skewsym(qmr) + 2*qmr*qmr')./4
end