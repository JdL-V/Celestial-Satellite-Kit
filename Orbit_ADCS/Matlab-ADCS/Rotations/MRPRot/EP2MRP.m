function var = EP2MRP(q)
    qmr = MRP(q.x(2:4)./(1 + q.x(1)), ...
                q.outFrame, q.inFrame);
    var = qmr;
end