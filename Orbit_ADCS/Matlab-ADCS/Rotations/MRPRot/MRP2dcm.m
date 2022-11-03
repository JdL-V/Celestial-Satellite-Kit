function var = MRP2dcm(qmr)
    qmrMat = skewsym(qmr.x);
    qmrSq = dot(qmr.x,qmr.x);
    dcm = DCM(eye(3, 3) + (8*qmrMat^2 - 4*(1 - qmrSq)*qmrMat)/(1 + qmrSq)^2, ...
                qmr.outFrame, qmr.inFrame);
    var = dcm;
end