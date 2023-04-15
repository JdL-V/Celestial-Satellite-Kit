function var = EP2CRP(q)
    if isstruct(q)
        vec = q.x;
        outFrame = q.outFrame;
        inFrame = q.inFrame;
    else
        checkVec(q, 4)
        vec = q;
        outFrame = "unk" + num2str(randi(5000));
        inFrame = "unk" + num2str(randi(5000));
    end
    
    qcr = CRP(vec(2:4)./vec(1), outFrame, inFrame);
    var = qcr;
end