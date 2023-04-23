function var = dcm2PRV(dcm)
    if isstruct(dcm)
        Mat = dcm.Mat;
        outFrame = dcm.outFrame;
        inFrame = dcm.inFrame;
    else
        checkMat(dcm, [3,3]);
        Mat = dcm;
        outFrame = "unk" + num2str(randi(5000));
        inFrame = "unk" + num2str(randi(5000));
    end

    Angle = acos(0.5*(Mat(1,1) + Mat(2,2) + Mat(3,3) - 1));
    x = normalize(1/(2*sin(Angle)).*[Mat(2,3)-Mat(3,2), Mat(3,1)-Mat(1,3), Mat(1,2)-Mat(2,1)], 'norm');
    prv = PRV(Angle, x, outFrame, inFrame);
    var = prv;
end