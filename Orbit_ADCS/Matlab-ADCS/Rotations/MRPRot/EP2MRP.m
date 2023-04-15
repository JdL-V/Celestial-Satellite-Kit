function var = EP2MRP(q)
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

    qmr = MRP(vec(2:4)./(1 + vec(1)), ...
                outFrame, inFrame);
    var = qmr;
end