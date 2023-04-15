function var = MRP2dcm(qmr)
    if isstruct(qmr)
        vec = qmr.x;
        outFrame = qmr.outFrame;
        inFrame = qmr.inFrame;
    else
        checkVec(q, 3)
        vec = qmr;
        outFrame = "unk" + num2str(randi(5000));
        inFrame = "unk" + num2str(randi(5000));
    end

    qmrMat = skewsym(vec);
    qmrSq = dot(vec,vec);
    dcm = DCM(eye(3, 3) + (8*qmrMat^2 - 4*(1 - qmrSq)*qmrMat)/(1 + qmrSq)^2, ...
                outFrame, inFrame);
    var = dcm;
end