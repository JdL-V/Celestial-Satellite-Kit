function var = CRP2EP(qcr)
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

    q = quat([1; vec(1:3)]./sqrt(1 + dot(vec,vec)), ...
                    outFrame, inFrame);
    var = q;
end