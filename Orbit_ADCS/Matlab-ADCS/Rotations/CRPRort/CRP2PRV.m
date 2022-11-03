function var = CRP2PRV(qcr)
    T = atan(qcr.x);
    e = normalize(T, 'norm');
    prv = PRV(2*T(1)/e(1), e, qcr.outFrame, qcr.inFrame);
    var = prv;
end