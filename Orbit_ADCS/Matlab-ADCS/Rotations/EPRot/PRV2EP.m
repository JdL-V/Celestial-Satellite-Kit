function var = PRV2EP(prv)
    q = quat([cos(prv.Angle/2), prv.x(1)*sin(prv.Angle/2), prv.x(2)*sin(prv.Angle/2), prv.x(3)*sin(prv.Angle/2)],...
                        prv.outFrame, prv.inFrame);
    var = q;
end