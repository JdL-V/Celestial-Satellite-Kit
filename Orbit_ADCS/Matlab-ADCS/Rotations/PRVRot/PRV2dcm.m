function var = PRV2dcm(prv)
    sig = 1 - cos(prv.Angle);
    dcm = DCM(zeros(3,3), prv.outFrame, prv.inFrame);
    for i = 1:3
        for j = 1:3
            dcm.Mat(i,j) = prv.x(i)*prv.x(j)*sig;
            if i == j
                dcm.Mat(i,j) = dcm.Mat(i,j) + cos(prv.Angle);
            else
                dcm.Mat(i,j) = dcm.Mat(i,j) + (-1)^(i + j + (i < j))*sin(prv.Angle)*prv.x([1 2 3]*prod([1 2 3]' ~= [i j], 2));
            end
        end
    end
    var = dcm;
end