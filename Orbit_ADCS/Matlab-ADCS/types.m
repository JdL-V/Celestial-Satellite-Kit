function var = types
    var.DCM         = @DCM;
    var.EulerAng    = @EulerAng;
    var.PRV         = @PRV;
    var.quaternion  = @quaternion;
    var.CRP         = @CRP;
    var.MRP         = @MRP;
end

function var = DCM(Mat, outFrame, inFrame)
    checkMat(Mat, [3,3])
    checkFrames(inFrame, outFrame)
    var = struct('Mat',      {Mat},      ...
                 'outFrame', {outFrame}, ...
                 'inFrame',  {inFrame});
end

function var = EulerAng(Angles, Sequence, outFrame, inFrame)
    checkVec(Angles, 3)
    checkVec(Sequence, 3)
    checkFrames(inFrame, outFrame)
    var = struct('Angles',   {Angles(:)},   ...
                 'Sequence', {Sequence(:)'}, ...
                 'outFrame', {outFrame}, ...
                 'inFrame',  {inFrame});
end

function var = PRV(Angle, Vector, outFrame, inFrame)
    checkVec(Angle, 1)
    checkVec(Vector, 3)
    checkFrames(inFrame, outFrame)
    var = struct('Angle',    {Angle},    ...
                 'x',        {Vector(:)},   ...
                 'outFrame', {outFrame}, ...
                 'inFrame',  {inFrame});
end

function var = quaternion(Vector, outFrame, inFrame)
    checkVec(Vector, 4)
    checkFrames(inFrame, outFrame)
    var = struct('x',        {Vector(:)},   ...
                 'outFrame', {outFrame}, ...
                 'inFrame',  {inFrame});
end

function var = CRP(Vector, outFrame, inFrame)
    checkVec(Vector, 3)
    checkFrames(inFrame, outFrame)
    var = struct('x',        {Vector(:)},   ...
                 'outFrame', {outFrame}, ...
                 'inFrame',  {inFrame});
end

function var = MRP(Vector, outFrame, inFrame)
    checkVec(Vector, 3)
    checkFrames(inFrame, outFrame)
    var = struct('x',        {Vector(:)},   ...
                 'outFrame', {outFrame}, ...
                 'inFrame',  {inFrame});
end

function var = FMatrix(Matrix, Frame)
    var = struct('Mat',      {Matrix},   ...
                 'Frame',    {Frame});
end

function var = FVector(Vector, Frame)
    var = struct('x',        {Vector(:)},   ...
                 'Frame',    {Frame});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function var = checkFrames(inFrame, outFrame)
    if ~(ischar(outFrame) || isstring(outFrame))
        error('Output frame is not a string')
    end

    if ~(ischar(inFrame) || isstring(inFrame))
        error('Input frame is not a string')
    end
end

function var = checkMat(Mat, Dim)
    if ~(prod(size(Mat) == Dim))
        error('Matrix dimesions not %d by %d', Dim(1), Dim(2))
    end
end

function var = checkVec(Vec, Dim)
    if ~(length(Vec(:)) == Dim)
        error('Vector dimesion not %d', Dim)
    end
end