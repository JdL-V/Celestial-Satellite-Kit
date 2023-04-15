function var = SumMRP(qmr1, qmr2)
    if qmr1.outFrame ~= qmr2.inFrame
        error("Rodrigues parameters are not frame-compatible")
    end

    qmr1Sq = dot(qmr1.x,qmr1.x);
    qmr2Sq = dot(qmr2.x,qmr2.x);
    qmr = MRP(((1 - qmr1Sq)*qmr2.x + (1 - qmr2Sq)*qmr1.x - 2*cross(qmr2.x,qmr1.x)) ...
                ./(1 + qmr1Sq*qmr2Sq - 2*dot(qmr1.x,qmr2.x)), ...
                qmr2.outFrame, qmr1.inFrame);
    var = qmr;
end