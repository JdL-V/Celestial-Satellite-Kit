function var = quat(Vector, outFrame, inFrame)
    checkVec(Vector, 4)
    checkFrames(inFrame, outFrame)
    var = struct('x',        {Vector(:)},   ...
                 'outFrame', {outFrame}, ...
                 'inFrame',  {inFrame});
end