function var = MRP2EP(qmr)
    sig = dot(qmr.x,qmr.x);
    q = quat([(1 - sig)/(1 + sig) 2*qmr.x(1:3)'/(1 + sig)], ...
                    qmr.outFrame, qmr.inFrame);
    var = q;
end