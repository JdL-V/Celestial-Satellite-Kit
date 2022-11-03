function var = dcm2PRV(dcm)
    Angle = acos(0.5*(dcm.Mat(1,1) + dcm.Mat(2,2) + dcm.Mat(3,3) - 1));
    x = 1/(2*sin(Angle)).*[dcm.Mat(2,3)-dcm.Mat(3,2), dcm.Mat(3,1)-dcm.Mat(1,3), dcm.Mat(1,2)-dcm.Mat(2,1)];
    prv = PRV(Angle, x, dcm.outFrame, dcm.inFrame);
    var = prv;
end