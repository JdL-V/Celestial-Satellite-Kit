function var = CRP2MRP(qcr)
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
    
    qmr = MRP(vec./(1 + sqrt(1 + dot(vec,vec))), ...
                outFrame, inFrame);
    var = qmr;
end