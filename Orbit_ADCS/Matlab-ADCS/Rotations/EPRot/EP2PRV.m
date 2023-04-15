function var = EP2PRV(q)
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

    angle = 2*acos(vec(1));
    prv = PRV(angle, vec(2:4)./sin(angle/2),...
                outFrame, inFrame);
    var = prv;
end