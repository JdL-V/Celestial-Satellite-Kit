% Converts direction cosine matrix to Euler angles
% dcm:     DCM type
% rot_seq: Integer Vector
function var = dcm2euler(dcm, rot_seq)
    % Switch between different rotation sequences
    % Generate Euler angles and save as an euler angle type
    if isstruct(dcm)
        Mat = dcm.Mat;
        outFrame = dcm.outFrame;
        inFrame = dcm.inFrame;
    else
        checkMat(dcm, [3,3])
        Mat = dcm;
        outFrame = "unk" + num2str(randi(5000));
        inFrame = "unk" + num2str(randi(5000));
    end

    switch rot_seq
    case 'ZYX'
        % Check for singularities.
        if ~(abs(Mat(1, 3)) >= 1 - eps())
            var = EulerAng([
                atanMk2(Mat(1, 2), Mat(1, 1)),
                asin(-Mat(1, 3)),
                atanMk2(Mat(2, 3), Mat(3, 3)),
                ], [3,2,1], outFrame, inFrame);
        else
            var = EulerAng([
                atanMk2(-Mat(2, 1), Mat(2, 2)),
                asinMk2(-Mat(1, 3)),
                0.,
                ], [3,2,1], outFrame, inFrame);
        end
    case 'XYX'
        % Check for singularities.
        if ~(abs(Mat(1, 1)) >= 1 - eps())
            var = EulerAng([
                atanMk2(Mat(1, 2), -Mat(1, 3)),
                acos(Mat(1, 1)),
                atanMk2(Mat(2, 1), Mat(3, 1)),
                ], [1,2,1], outFrame, inFrame);
        else
            var = EulerAng([
                atanMk2(Mat(2, 3), Mat(2, 2)),
                acosMk2(Mat(1, 1)),
                0.,
                ], [1,2,1], outFrame, inFrame);
        end
    case 'XYZ'
        % Check for singularities.
        if ~(abs(Mat(3, 1)) >= 1 - eps())
            var = EulerAng([
                atanMk2(-Mat(3, 2), Mat(3, 3)),
                asin(Mat(3, 1)),
                atanMk2(-Mat(2, 1), Mat(1, 1)),
                ], [1,2,3], outFrame, inFrame);
        else
            var = EulerAng([
                atanMk2(Mat(2, 3), Mat(2, 2)),
                asinMk2(Mat(3, 1)),
                0.,
                ], [1,2,3], outFrame, inFrame);
        end
    case 'XZX'
        % Check for singularities.
        if ~(abs(Mat(1, 1)) >= 1 - eps())
            var = EulerAng([
                atanMk2(Mat(1, 3), Mat(1, 2)),
                acos(Mat(1, 1)),
                atanMk2(Mat(3, 1), -Mat(2, 1)),
                ], [1,2,3], outFrame, inFrame);
        else
            var = EulerAng([
                atanMk2(-Mat(3, 2), Mat(3, 3)),
                acosMk2(Mat(1, 1)),
                0.,
                ], [1,2,3], outFrame, inFrame);
        end
    case 'XZY'
        % Check for singularities.
        if ~(abs(Mat(2, 1)) >= 1 - eps())
            var = EulerAng([
                atanMk2(Mat(2, 3), Mat(2, 2)),
                asin(-Mat(2, 1)),
                atanMk2(Mat(3, 1), Mat(1, 1)),
                ], [1,3,2], outFrame, inFrame);
        else
            var = EulerAng([
                atanMk2(-Mat(3, 2), Mat(3, 3)),
                asinMk2(-Mat(2, 1)),
                0.,
                ], [1,3,2], outFrame, inFrame);
        end
    case 'YXY'
        % Check for singularities.
        if ~(abs(Mat(2, 2)) >= 1 - eps())
            var = EulerAng([
                atanMk2(Mat(2, 1), Mat(2, 3)),
                acos(Mat(2, 2)),
                atanMk2(Mat(1, 2), -Mat(3, 2)),
                ], [2,1,2], outFrame, inFrame);
        else
            var = EulerAng([
                atanMk2(-Mat(1, 3), Mat(1, 1)),
                acosMk2(Mat(2, 2)),
                0.,
                ], [2,1,2], outFrame, inFrame);
        end
    case 'YXZ'
        if ~(abs(Mat(3, 2)) >= 1 - eps())
            var = EulerAng([
                atanMk2(Mat(3, 1), Mat(3, 3)),
                asin(-Mat(3, 2)),
                atanMk2(Mat(1, 2), Mat(2, 2)),
                ], [2,1,3], outFrame, inFrame);
        else
            var = EulerAng([
                atanMk2(-Mat(1, 3), Mat(1, 1)),
                asinMk2(-Mat(3, 2)),
                0.,
                ], [2,1,3], outFrame, inFrame);
        end
    case 'YZX'
        % Check for singularities.
        if ~(abs(Mat(1, 2)) >= 1 - eps())
            var = EulerAng([
                atanMk2(-Mat(1, 3), Mat(1, 1)),
                asin(Mat(1, 2)),
                atanMk2(-Mat(3, 2), Mat(2, 2)),
                ], [2,3,1], outFrame, inFrame);
        else
            var = EulerAng([
                atanMk2(Mat(3, 1), Mat(3, 3)),
                asinMk2(Mat(1, 2)),
                0.,
                ], [2,3,1], outFrame, inFrame);
        end
    case 'YZY'
        % Check for singularities.
        if ~(abs(Mat(2, 2)) >= 1 - eps())
            var = EulerAng([
                atanMk2(Mat(2, 3), -Mat(2, 1)),
                acos(Mat(2, 2)),
                atanMk2(Mat(3, 2), Mat(1, 2)),
                ], [2,3,2], outFrame, inFrame);
        else
            var = EulerAng([
                atanMk2(Mat(3, 1), Mat(3, 3)),
                acosMk2(Mat(2, 2)),
                0.,
                ], [2,3,2], outFrame, inFrame);
        end
    case 'ZXY'
        % Check for singularities.
        if ~(abs(Mat(2, 3)) >= 1 - eps())
            var = EulerAng([
                atanMk2(-Mat(2, 1), Mat(2, 2)),
                asin(Mat(2, 3)),
                atanMk2(-Mat(1, 3), Mat(3, 3)),
                ], [3,1,2], outFrame, inFrame);
        else
            var = EulerAng([
                atanMk2(Mat(1, 2), Mat(1, 1)),
                asinMk2(Mat(2, 3)),
                0.,
                ], [3,1,2], outFrame, inFrame);
        end
    case 'ZXZ'
        % Check for singularities.
        if ~(abs(Mat(3, 3)) >= 1 - eps())
            var = EulerAng([
                atanMk2(Mat(3, 1), -Mat(3, 2)),
                acos(Mat(3, 3)),
                atanMk2(Mat(1, 3), Mat(2, 3)),
                ], [3,1,3], outFrame, inFrame);
        else
            var = EulerAng([
                atanMk2(Mat(1, 2), Mat(1, 1)),
                acosMk2(Mat(3, 3)),
                0.,
                ], [3,1,3], outFrame, inFrame);
        end
    case 'ZYZ'
        % Check for singularities.
        if ~(abs(Mat(3, 3)) >= 1 - eps())
            var = EulerAng([
                atanMk2(Mat(3, 2), Mat(3, 1)),
                acos(Mat(3, 3)),
                atanMk2(Mat(2, 3), -Mat(1, 3)),
                ], [3,2,3], outFrame, inFrame);
        else
            var = EulerAng([
                atanMk2(-Mat(2, 1), Mat(2, 2)),
                acosMk2(Mat(3, 3)),
                0.,
                ], [3,2,3], outFrame, inFrame);
        end
    otherwise
        error("The rotation sequence is not valid.")
    end
end