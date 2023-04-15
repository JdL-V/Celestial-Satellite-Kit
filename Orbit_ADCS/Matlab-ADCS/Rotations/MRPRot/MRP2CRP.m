function var = MRP2CRP(qmr)
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

    qcr = CRP(vec./(1 - dot(vec,vec)), ...
                outFrame, inFrame);
    var = qcr;
end