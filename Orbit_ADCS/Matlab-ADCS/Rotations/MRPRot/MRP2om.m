function var = MRP2om(qmr)
    var = (4/(1 + dot(qmr,qmr))^2)*(om2MRP(qmr)');
end