function var = PRV2MRP(prv)
    qmr = MRP(tan(prv.Angle/4).*prv.x, ...
                prv.outFrame, prv.inFrame);
    var = qmr;
end