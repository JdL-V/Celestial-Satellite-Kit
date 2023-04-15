function var = PRV2MRP(prv)
    if isstruct(prv)
        vec = prv.x;
        angle = prv.Angle
        outFrame = prv.outFrame;
        inFrame = prv.inFrame;
    else
        checkVec(prv, 4)
        angle = prv(1)
        vec = prv(2:4);
        outFrame = "unk" + num2str(randi(5000));
        inFrame = "unk" + num2str(randi(5000));
    end

    qmr = MRP(tan(angle/4).*vec, ...
                outFrame, inFrame);
    var = qmr;
end