function var = DCMtrs(dcm)
    var = DCM(dcm.Mat', dcm.inFrame, dcm.outFrame);
end