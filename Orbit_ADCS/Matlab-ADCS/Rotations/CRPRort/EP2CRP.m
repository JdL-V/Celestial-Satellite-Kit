function var = EP2CRP(q)
    qcr = CRP(q.x(2:4)./q.x(1), q.outFrame, q.inFrame);
    var = qcr;
end