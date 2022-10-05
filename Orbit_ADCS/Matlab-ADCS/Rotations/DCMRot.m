function var = DCMRot
    global TP 
    TP = types;

    var.DCMmul    = @DCMmul;
    var.DCMtrs    = @DCMtrs
end

function var = DCMmul(dcm2, dcm1)
    if dcm1.outFrame != dcm2.inFrame
        error("dcms are not frame-compatible")
    end
    var = TP.DCM(dcm2.Mat*dcm1.Mat, dcm2.outFrame, dcm1.inFrame)
end

function var = DCMtrs(dcm)
    var = TP.DCM(dcm.Mat', dcm.inFrame, dcm.outFrame);
end