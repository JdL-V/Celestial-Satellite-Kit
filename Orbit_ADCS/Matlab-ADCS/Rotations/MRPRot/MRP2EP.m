function var = MRP2EP(qmr)
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

    sig = dot(vec, vec);
    q = quat([(1 - sig)/(1 + sig) 2*vec(1:3)'/(1 + sig)], ...
                    outFrame, inFrame);
    var = q;
end