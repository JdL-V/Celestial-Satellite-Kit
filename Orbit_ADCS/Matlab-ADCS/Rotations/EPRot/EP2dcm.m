function var = EP2dcm(q)
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

    dcm = DCM([vec(1)^2 + vec(2)^2 - vec(3)^2 - vec(4)^2   2*(vec(2)*vec(3) + vec(1)*vec(4))           2*(vec(2)*vec(4) - vec(1)*vec(3))
               2*(vec(2)*vec(3) - vec(1)*vec(4))           vec(1)^2 - vec(2)^2 + vec(3)^2 - vec(4)^2   2*(vec(3)*vec(4) + vec(1)*vec(2)) 
               2*(vec(2)*vec(4) + vec(1)*vec(3))           2*(vec(3)*vec(4) - vec(1)*vec(2))           vec(1)^2 - vec(2)^2 - vec(3)^2 + vec(4)^2],...
               outFrame, inFrame);
    var = dcm;
end