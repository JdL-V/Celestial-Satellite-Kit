function var = DCMmul(dcm2, dcm1)
    if dcm1.outFrame ~= dcm2.inFrame
        error("dcms are not frame-compatible")
    end
    var = DCM(dcm2.Mat*dcm1.Mat, dcm2.outFrame, dcm1.inFrame);
end