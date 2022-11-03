function var = CRP2EP(qcr)
    q = quat([1, qcr.x(1:3)]./sqrt(1 + dot(qcr,qcr)), ...
                    qcr.outFrame, qcr.inFrame);
    var = q;
end