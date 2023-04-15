function var = PRV2dcm(prv)
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

    sig = 1 - cos(angle);
    dcm = DCM(zeros(3,3), outFrame, inFrame);
    for i = 1:3
        for j = 1:3
            dcm.Mat(i,j) = vec(i)*vec(j)*sig;
            if i == j
                dcm.Mat(i,j) = dcm.Mat(i,j) + cos(angle);
            else
                dcm.Mat(i,j) = dcm.Mat(i,j) + (-1)^(i + j + (i < j))*sin(angle)*vec([1 2 3]*prod([1 2 3]' ~= [i j], 2));
            end
        end
    end
    var = dcm;
end