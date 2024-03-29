function var = SumCRP(qcr1, qcr2)
    if qcr1.outFrame ~= qcr2.inFrame
        error("Rodrigues parameters are not frame-compatible")
    end

    qcr = CRP((qcr2.x + qcr1.x - cross(qcr2.x,qcr1.x))./(1 - dot(qcr2.x,qcr1.x)), ...
                qcr2.outFrame, qcr1.inFrame);
    var = qcr;
end