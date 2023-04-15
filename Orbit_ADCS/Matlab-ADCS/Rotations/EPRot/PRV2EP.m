function var = PRV2EP(prv)
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

    q = quat([cos(angle/2), vec(1)*sin(angle/2), vec(2)*sin(angle/2), vec(3)*sin(angle/2)],...
                        outFrame, inFrame);
    var = q;
end