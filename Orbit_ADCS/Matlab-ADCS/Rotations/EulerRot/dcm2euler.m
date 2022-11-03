% Converts direction cosine matrix to Euler angles
% dcm:     DCM type
% rot_seq: Integer Vector
function var = dcm2euler(dcm, rot_seq)
    % Switch between different rotation sequences
    % Generate Euler angles and save as an euler angle type
    switch rot_seq
    case 'ZYX'
        % Check for singularities.
        if ~(abs(dcm.Mat(1, 3)) >= 1 - eps())
            var = EulerAng([
                atanMk2(dcm.Mat(1, 2), dcm.Mat(1, 1)),
                asin(-dcm.Mat(1, 3)),
                atanMk2(dcm.Mat(2, 3), dcm.Mat(3, 3)),
                ], [3,2,1], dcm.outFrame, dcm.inFrame);
        else
            var = EulerAng([
                atanMk2(-dcm.Mat(2, 1), dcm.Mat(2, 2)),
                asinMk2(-dcm.Mat(1, 3)),
                0.,
                ], [3,2,1], dcm.outFrame, dcm.inFrame);
        end
    case 'XYX'
        % Check for singularities.
        if ~(abs(dcm.Mat(1, 1)) >= 1 - eps())
            var = EulerAng([
                atanMk2(dcm.Mat(1, 2), -dcm.Mat(1, 3)),
                acos(dcm.Mat(1, 1)),
                atanMk2(dcm.Mat(2, 1), dcm.Mat(3, 1)),
                ], [1,2,1], dcm.outFrame, dcm.inFrame);
        else
            var = EulerAng([
                atanMk2(dcm.Mat(2, 3), dcm.Mat(2, 2)),
                acosMk2(dcm.Mat(1, 1)),
                0.,
                ], [1,2,1], dcm.outFrame, dcm.inFrame);
        end
    case 'XYZ'
        % Check for singularities.
        if ~(abs(dcm.Mat(3, 1)) >= 1 - eps())
            var = EulerAng([
                atanMk2(-dcm.Mat(3, 2), dcm.Mat(3, 3)),
                asin(dcm.Mat(3, 1)),
                atanMk2(-dcm.Mat(2, 1), dcm.Mat(1, 1)),
                ], [1,2,3], dcm.outFrame, dcm.inFrame);
        else
            var = EulerAng([
                atanMk2(dcm.Mat(2, 3), dcm.Mat(2, 2)),
                asinMk2(dcm.Mat(3, 1)),
                0.,
                ], [1,2,3], dcm.outFrame, dcm.inFrame);
        end
    case 'XZX'
        % Check for singularities.
        if ~(abs(dcm.Mat(1, 1)) >= 1 - eps())
            var = EulerAng([
                atanMk2(dcm.Mat(1, 3), dcm.Mat(1, 2)),
                acos(dcm.Mat(1, 1)),
                atanMk2(dcm.Mat(3, 1), -dcm.Mat(2, 1)),
                ], [1,2,3], dcm.outFrame, dcm.inFrame);
        else
            var = EulerAng([
                atanMk2(-dcm.Mat(3, 2), dcm.Mat(3, 3)),
                acosMk2(dcm.Mat(1, 1)),
                0.,
                ], [1,2,3], dcm.outFrame, dcm.inFrame);
        end
    case 'XZY'
        % Check for singularities.
        if ~(abs(dcm.Mat(2, 1)) >= 1 - eps())
            var = EulerAng([
                atanMk2(dcm.Mat(2, 3), dcm.Mat(2, 2)),
                asin(-dcm.Mat(2, 1)),
                atanMk2(dcm.Mat(3, 1), dcm.Mat(1, 1)),
                ], [1,3,2], dcm.outFrame, dcm.inFrame);
        else
            var = EulerAng([
                atanMk2(-dcm.Mat(3, 2), dcm.Mat(3, 3)),
                asinMk2(-dcm.Mat(2, 1)),
                0.,
                ], [1,3,2], dcm.outFrame, dcm.inFrame);
        end
    case 'YXY'
        % Check for singularities.
        if ~(abs(dcm.Mat(2, 2)) >= 1 - eps())
            var = EulerAng([
                atanMk2(dcm.Mat(2, 1), dcm.Mat(2, 3)),
                acos(dcm.Mat(2, 2)),
                atanMk2(dcm.Mat(1, 2), -dcm.Mat(3, 2)),
                ], [2,1,2], dcm.outFrame, dcm.inFrame);
        else
            var = EulerAng([
                atanMk2(-dcm.Mat(1, 3), dcm.Mat(1, 1)),
                acosMk2(dcm.Mat(2, 2)),
                0.,
                ], [2,1,2], dcm.outFrame, dcm.inFrame);
        end
    case 'YXZ'
        if ~(abs(dcm.Mat(3, 2)) >= 1 - eps())
            var = EulerAng([
                atanMk2(dcm.Mat(3, 1), dcm.Mat(3, 3)),
                asin(-dcm.Mat(3, 2)),
                atanMk2(dcm.Mat(1, 2), dcm.Mat(2, 2)),
                ], [2,1,3], dcm.outFrame, dcm.inFrame);
        else
            var = EulerAng([
                atanMk2(-dcm.Mat(1, 3), dcm.Mat(1, 1)),
                asinMk2(-dcm.Mat(3, 2)),
                0.,
                ], [2,1,3], dcm.outFrame, dcm.inFrame);
        end
    case 'YZX'
        % Check for singularities.
        if ~(abs(dcm.Mat(1, 2)) >= 1 - eps())
            var = EulerAng([
                atanMk2(-dcm.Mat(1, 3), dcm.Mat(1, 1)),
                asin(dcm.Mat(1, 2)),
                atanMk2(-dcm.Mat(3, 2), dcm.Mat(2, 2)),
                ], [2,3,1], dcm.outFrame, dcm.inFrame);
        else
            var = EulerAng([
                atanMk2(dcm.Mat(3, 1), dcm.Mat(3, 3)),
                asinMk2(dcm.Mat(1, 2)),
                0.,
                ], [2,3,1], dcm.outFrame, dcm.inFrame);
        end
    case 'YZY'
        % Check for singularities.
        if ~(abs(dcm.Mat(2, 2)) >= 1 - eps())
            var = EulerAng([
                atanMk2(dcm.Mat(2, 3), -dcm.Mat(2, 1)),
                acos(dcm.Mat(2, 2)),
                atanMk2(dcm.Mat(3, 2), dcm.Mat(1, 2)),
                ], [2,3,2], dcm.outFrame, dcm.inFrame);
        else
            var = EulerAng([
                atanMk2(dcm.Mat(3, 1), dcm.Mat(3, 3)),
                acosMk2(dcm.Mat(2, 2)),
                0.,
                ], [2,3,2], dcm.outFrame, dcm.inFrame);
        end
    case 'ZXY'
        % Check for singularities.
        if ~(abs(dcm.Mat(2, 3)) >= 1 - eps())
            var = EulerAng([
                atanMk2(-dcm.Mat(2, 1), dcm.Mat(2, 2)),
                asin(dcm.Mat(2, 3)),
                atanMk2(-dcm.Mat(1, 3), dcm.Mat(3, 3)),
                ], [3,1,2], dcm.outFrame, dcm.inFrame);
        else
            var = EulerAng([
                atanMk2(dcm.Mat(1, 2), dcm.Mat(1, 1)),
                asinMk2(dcm.Mat(2, 3)),
                0.,
                ], [3,1,2], dcm.outFrame, dcm.inFrame);
        end
    case 'ZXZ'
        % Check for singularities.
        if ~(abs(dcm.Mat(3, 3)) >= 1 - eps())
            var = EulerAng([
                atanMk2(dcm.Mat(3, 1), -dcm.Mat(3, 2)),
                acos(dcm.Mat(3, 3)),
                atanMk2(dcm.Mat(1, 3), dcm.Mat(2, 3)),
                ], [3,1,3], dcm.outFrame, dcm.inFrame);
        else
            var = EulerAng([
                atanMk2(dcm.Mat(1, 2), dcm.Mat(1, 1)),
                acosMk2(dcm.Mat(3, 3)),
                0.,
                ], [3,1,3], dcm.outFrame, dcm.inFrame);
        end
    case 'ZYZ'
        % Check for singularities.
        if ~(abs(dcm.Mat(3, 3)) >= 1 - eps())
            var = EulerAng([
                atanMk2(dcm.Mat(3, 2), dcm.Mat(3, 1)),
                acos(dcm.Mat(3, 3)),
                atanMk2(dcm.Mat(2, 3), -dcm.Mat(1, 3)),
                ], [3,2,3], dcm.outFrame, dcm.inFrame);
        else
            var = EulerAng([
                atanMk2(-dcm.Mat(2, 1), dcm.Mat(2, 2)),
                acosMk2(dcm.Mat(3, 3)),
                0.,
                ], [3,2,3], dcm.outFrame, dcm.inFrame);
        end
    otherwise
        error("The rotation sequence is not valid.")
    end
end