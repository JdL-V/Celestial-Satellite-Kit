function var = MRP2PRV(qmr)
    T = atan(qmr.x);
    e = normalize(T, 'norm');
    prv = PRV(4*T(1)/e(1), e, qmr.outFrame, qmr.inFrame);
    var = prv;
end