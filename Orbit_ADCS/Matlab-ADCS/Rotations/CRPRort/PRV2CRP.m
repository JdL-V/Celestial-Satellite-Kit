function var = PRV2CRP(prv)
    qcr = CRP(tan(prv.Angle/2).*prv.x, prv.outFrame, prv.inFrame);
    var = qcr;
end