function var = MRP2PRV(qmr)
    if isstruct(qmr)
        vec = qmr.x;
        outFrame = qmr.outFrame;
        inFrame = qmr.inFrame;
    else
        checkVec(q, 3)
        vec = qmr;
        outFrame = "unk" + num2str(randi(5000));
        inFrame = "unk" + num2str(randi(5000));
    end
    
    T = atan(vec);
    e = normalize(T, 'norm');
    prv = PRV(4*T(1)/e(1), e, outFrame, inFrame);
    var = prv;
end