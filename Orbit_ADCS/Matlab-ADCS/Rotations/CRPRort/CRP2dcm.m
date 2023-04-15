function var = CRP2dcm(qcr)
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

    dcm = DCM([ 1 + vec(1)^2 - vec(2)^2 - vec(3)^2    2*(vec(1)*vec(2) + vec(3))            2*(vec(1)*vec(3) - vec(2))
                2*(vec(1)*vec(2) - vec(3))            1 - vec(1)^2 + vec(2)^2 - vec(3)^2    2*(vec(2)*vec(3) + vec(1)) 
                2*(vec(1)*vec(3) + vec(2))            2*(vec(2)*vec(3) - vec(1))            1 - vec(1)^2 - vec(2)^2 + vec(3)^2]./(1 + dot(vec,vec)), ...
            outFrame, inFrame);
    var = dcm;
    % var = ((1 - dot(qcr,qcr))*Matrix(I, 3, 3) + 2*qcr*qcr' - 2*skewsym(qcr))./(1 + dot(qcr,qcr))
end