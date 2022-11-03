function var = DCM(Mat, outFrame, inFrame)
    checkMat(Mat, [3,3])
    checkFrames(inFrame, outFrame)
    var = struct('Mat',      {Mat},      ...
                 'outFrame', {outFrame}, ...
                 'inFrame',  {inFrame});
end