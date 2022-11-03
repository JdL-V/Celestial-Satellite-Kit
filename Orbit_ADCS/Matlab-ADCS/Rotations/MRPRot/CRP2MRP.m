function var = CRP2MRP(qcr)
    qmr = MRP(qcr.x./(1 + sqrt(1 + dot(qcr.x,qcr.x))), ...
                qcr.outFrame, qcr.inFrame);
    var = qmr;
end