function var = CRP2om(qcr)
    var = 2/(1 + dot(qcr, qcr)).*(eye(3,3) - skewsym(qcr));
end