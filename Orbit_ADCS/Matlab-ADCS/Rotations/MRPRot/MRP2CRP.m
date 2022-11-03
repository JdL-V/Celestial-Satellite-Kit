function var = MRP2CRP(qmr)
    qcr = CRP(qmr.x./(1 - dot(qmr.x,qmr.x)), ...
                qmr.outFrame, qmr.inFrame);
    var = qcr;
end