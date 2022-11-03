function var = MRP(Vector, outFrame, inFrame)
    checkVec(Vector, 3)
    checkFrames(inFrame, outFrame)
    var = struct('x',        {Vector(:)},   ...
                 'outFrame', {outFrame}, ...
                 'inFrame',  {inFrame});
end