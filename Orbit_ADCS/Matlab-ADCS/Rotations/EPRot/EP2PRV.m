function var = EP2PRV(q)
    angle = 2*acos(q.x(1));
    prv = PRV(angle, q.x(2:4)./sin(angle/2),...
                q.outFrame, q.inFrame);
    var = prv;
end