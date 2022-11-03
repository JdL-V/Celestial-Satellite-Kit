function var = SumEP(q1, q2)
    q = quat([q2.x(1) -q2.x(2) -q2.x(3) -q2.x(4);
              q2.x(2)  q2.x(1)  q2.x(4) -q2.x(3);
              q2.x(3) -q2.x(4)  q2.x(1)  q2.x(2);
              q2.x(4)  q2.x(3) -q2.x(2)  q2.x(1)]*(q1.x),...
              q2.outFrame, q1.inFrame);
    var = q;
end