function var = CRP2PRV(qcr)
    if isstruct(qcr)
        vec = qcr.x;
        outFrame = qcr.outFrame;
        inFrame = qcr.inFrame;
    else
        checkVec(q, 3)
        vec = qcr;
        outFrame = "unk" + num2str(randi(5000));
        inFrame = "unk" + num2str(randi(5000));
    end
    
    T = atan(vec);
    e = normalize(T, 'norm');
    prv = PRV(2*T(1)/e(1), e, outFrame, inFrame);
    var = prv;
end