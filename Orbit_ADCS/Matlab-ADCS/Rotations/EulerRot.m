% Declaration of Function containing exportable functions
function var = EulerRot
    global TP 
    TP = types;

    var.euler2dcm   = @euler2dcm;
    var.dcm2euler   = @dcm2euler;
    var.euler2om    = @euler2om;
    var.om2euler    = @om2euler;
end

% Converts Euler angles to direction cosine matrix
% EA: EA type
function var = euler2dcm(EA)
    global TP
    Mat = eye(3,3);
    for ind = 1:length(EA.Sequence)
        % Multiply the rotation matrices
        Mat = RotMat(EA.Angles(ind), EA.Sequence(ind))*Mat;
    end
    % Save the resulting matrix as a DCM type
    var = TP.DCM(Mat, EA.outFrame, EA.inFrame);
end

% Generate rotation matrix by choosing the axis and angle of rotation
% Angle: Float
% Axis:  Integer
function var = RotMat(Angle, Axis) 
    if Axis == 1
        % First axis rotation
        var = [1   0           0
               0   cos(Angle)  sin(Angle)
               0  -sin(Angle)  cos(Angle)];
    elseif Axis == 2
        % Second axis rotation
        var = [cos(Angle)  0  -sin(Angle)
               0           1   0
               sin(Angle)  0   cos(Angle)];
    elseif Axis == 3
        % Third axis rotation
        var = [cos(Angle)  sin(Angle)  0
              -sin(Angle)  cos(Angle)  0
               0           0           1];
    end
end

% Converts direction cosine matrix to Euler angles
% dcm:     DCM type
% rot_seq: Integer Vector
function var = dcm2euler(dcm, rot_seq)
    global TP
    % Switch between different rotation sequences
    % Generate Euler angles and save as an euler angle type
    switch rot_seq
    case 'ZYX'
        % Check for singularities.
        if ~(abs(dcm.Mat(1, 3)) >= 1 - eps())
            var = TP.EulerAng([
                atanMk2(dcm.Mat(1, 2), dcm.Mat(1, 1)),
                asin(-dcm.Mat(1, 3)),
                atanMk2(dcm.Mat(2, 3), dcm.Mat(3, 3)),
                ], [3,2,1], dcm.outFrame, dcm.inFrame);
        else
            var = TP.EulerAng([
                atanMk2(-dcm.Mat(2, 1), dcm.Mat(2, 2)),
                asinMk2(-dcm.Mat(1, 3)),
                0.,
                ], [3,2,1], dcm.outFrame, dcm.inFrame);
        end
    case 'XYX'
        % Check for singularities.
        if ~(abs(dcm.Mat(1, 1)) >= 1 - eps())
            var = TP.EulerAng([
                atanMk2(dcm.Mat(1, 2), -dcm.Mat(1, 3)),
                acos(dcm.Mat(1, 1)),
                atanMk2(dcm.Mat(2, 1), dcm.Mat(3, 1)),
                ], [1,2,1], dcm.outFrame, dcm.inFrame);
        else
            var = TP.EulerAng([
                atanMk2(dcm.Mat(2, 3), dcm.Mat(2, 2)),
                acosMk2(dcm.Mat(1, 1)),
                0.,
                ], [1,2,1], dcm.outFrame, dcm.inFrame);
        end
    case 'XYZ'
        % Check for singularities.
        if ~(abs(dcm.Mat(3, 1)) >= 1 - eps())
            var = TP.EulerAng([
                atanMk2(-dcm.Mat(3, 2), dcm.Mat(3, 3)),
                asin(dcm.Mat(3, 1)),
                atanMk2(-dcm.Mat(2, 1), dcm.Mat(1, 1)),
                ], [1,2,3], dcm.outFrame, dcm.inFrame);
        else
            var = TP.EulerAng([
                atanMk2(dcm.Mat(2, 3), dcm.Mat(2, 2)),
                asinMk2(dcm.Mat(3, 1)),
                0.,
                ], [1,2,3], dcm.outFrame, dcm.inFrame);
        end
    case 'XZX'
        % Check for singularities.
        if ~(abs(dcm.Mat(1, 1)) >= 1 - eps())
            var = TP.EulerAng([
                atanMk2(dcm.Mat(1, 3), dcm.Mat(1, 2)),
                acos(dcm.Mat(1, 1)),
                atanMk2(dcm.Mat(3, 1), -dcm.Mat(2, 1)),
                ], [1,2,3], dcm.outFrame, dcm.inFrame);
        else
            var = TP.EulerAng([
                atanMk2(-dcm.Mat(3, 2), dcm.Mat(3, 3)),
                acosMk2(dcm.Mat(1, 1)),
                0.,
                ], [1,2,3], dcm.outFrame, dcm.inFrame);
        end
    case 'XZY'
        % Check for singularities.
        if ~(abs(dcm.Mat(2, 1)) >= 1 - eps())
            var = TP.EulerAng([
                atanMk2(dcm.Mat(2, 3), dcm.Mat(2, 2)),
                asin(-dcm.Mat(2, 1)),
                atanMk2(dcm.Mat(3, 1), dcm.Mat(1, 1)),
                ], [1,3,2], dcm.outFrame, dcm.inFrame);
        else
            var = TP.EulerAng([
                atanMk2(-dcm.Mat(3, 2), dcm.Mat(3, 3)),
                asinMk2(-dcm.Mat(2, 1)),
                0.,
                ], [1,3,2], dcm.outFrame, dcm.inFrame);
        end
    case 'YXY'
        % Check for singularities.
        if ~(abs(dcm.Mat(2, 2)) >= 1 - eps())
            var = TP.EulerAng([
                atanMk2(dcm.Mat(2, 1), dcm.Mat(2, 3)),
                acos(dcm.Mat(2, 2)),
                atanMk2(dcm.Mat(1, 2), -dcm.Mat(3, 2)),
                ], [2,1,2], dcm.outFrame, dcm.inFrame);
        else
            var = TP.EulerAng([
                atanMk2(-dcm.Mat(1, 3), dcm.Mat(1, 1)),
                acosMk2(dcm.Mat(2, 2)),
                0.,
                ], [2,1,2], dcm.outFrame, dcm.inFrame);
        end
    case 'YXZ'
        if ~(abs(dcm.Mat(3, 2)) >= 1 - eps())
            var = TP.EulerAng([
                atanMk2(dcm.Mat(3, 1), dcm.Mat(3, 3)),
                asin(-dcm.Mat(3, 2)),
                atanMk2(dcm.Mat(1, 2), dcm.Mat(2, 2)),
                ], [2,1,3], dcm.outFrame, dcm.inFrame);
        else
            var = TP.EulerAng([
                atanMk2(-dcm.Mat(1, 3), dcm.Mat(1, 1)),
                asinMk2(-dcm.Mat(3, 2)),
                0.,
                ], [2,1,3], dcm.outFrame, dcm.inFrame);
        end
    case 'YZX'
        % Check for singularities.
        if ~(abs(dcm.Mat(1, 2)) >= 1 - eps())
            var = TP.EulerAng([
                atanMk2(-dcm.Mat(1, 3), dcm.Mat(1, 1)),
                asin(dcm.Mat(1, 2)),
                atanMk2(-dcm.Mat(3, 2), dcm.Mat(2, 2)),
                ], [2,3,1], dcm.outFrame, dcm.inFrame);
        else
            var = TP.EulerAng([
                atanMk2(dcm.Mat(3, 1), dcm.Mat(3, 3)),
                asinMk2(dcm.Mat(1, 2)),
                0.,
                ], [2,3,1], dcm.outFrame, dcm.inFrame);
        end
    case 'YZY'
        % Check for singularities.
        if ~(abs(dcm.Mat(2, 2)) >= 1 - eps())
            var = TP.EulerAng([
                atanMk2(dcm.Mat(2, 3), -dcm.Mat(2, 1)),
                acos(dcm.Mat(2, 2)),
                atanMk2(dcm.Mat(3, 2), dcm.Mat(1, 2)),
                ], [2,3,2], dcm.outFrame, dcm.inFrame);
        else
            var = TP.EulerAng([
                atanMk2(dcm.Mat(3, 1), dcm.Mat(3, 3)),
                acosMk2(dcm.Mat(2, 2)),
                0.,
                ], [2,3,2], dcm.outFrame, dcm.inFrame);
        end
    case 'ZXY'
        % Check for singularities.
        if ~(abs(dcm.Mat(2, 3)) >= 1 - eps())
            var = TP.EulerAng([
                atanMk2(-dcm.Mat(2, 1), dcm.Mat(2, 2)),
                asin(dcm.Mat(2, 3)),
                atanMk2(-dcm.Mat(1, 3), dcm.Mat(3, 3)),
                ], [3,1,2], dcm.outFrame, dcm.inFrame);
        else
            var = TP.EulerAng([
                atanMk2(dcm.Mat(1, 2), dcm.Mat(1, 1)),
                asinMk2(dcm.Mat(2, 3)),
                0.,
                ], [3,1,2], dcm.outFrame, dcm.inFrame);
        end
    case 'ZXZ'
        % Check for singularities.
        if ~(abs(dcm.Mat(3, 3)) >= 1 - eps())
            var = TP.EulerAng([
                atanMk2(dcm.Mat(3, 1), -dcm.Mat(3, 2)),
                acos(dcm.Mat(3, 3)),
                atanMk2(dcm.Mat(1, 3), dcm.Mat(2, 3)),
                ], [3,1,3], dcm.outFrame, dcm.inFrame);
        else
            var = TP.EulerAng([
                atanMk2(dcm.Mat(1, 2), dcm.Mat(1, 1)),
                acosMk2(dcm.Mat(3, 3)),
                0.,
                ], [3,1,3], dcm.outFrame, dcm.inFrame);
        end
    case 'ZYZ'
        % Check for singularities.
        if ~(abs(dcm.Mat(3, 3)) >= 1 - eps())
            var = TP.EulerAng([
                atanMk2(dcm.Mat(3, 2), dcm.Mat(3, 1)),
                acos(dcm.Mat(3, 3)),
                atanMk2(dcm.Mat(2, 3), -dcm.Mat(1, 3)),
                ], [3,2,3], dcm.outFrame, dcm.inFrame);
        else
            var = TP.EulerAng([
                atanMk2(-dcm.Mat(2, 1), dcm.Mat(2, 2)),
                acosMk2(dcm.Mat(3, 3)),
                0.,
                ], [3,2,3], dcm.outFrame, dcm.inFrame);
        end
    otherwise
        error("The rotation sequence is not valid.")
    end
end

% Computes acos(x)        if |x| <= 1
% Computes acos(sign(x))  if |x| > 1
function var = acosMk2(x)
    if x > 1
        var = 0.;
    elseif x < -1
        var = pi;
    else
        var = acos(x);
    end
end

% Computes asin(x)        if |x| <= 1
% Computes asin(sign(x))  if |x| > 1
function var = asinMk2(x) 
    if x > 1
        var = +(pi / 2);
    elseif x < -1
        var = -(pi / 2);
    else
        var = asin(x);
    end
end

% Computes the Euler kinematic equation to transform 
% Euler Angle rates to Body rotation rates
% EA: EA type
function var = euler2om(EA)

    % Compute the sines and cosines.
    [s2, c2] = sincos(EA.Angles(2));
    [s3, c3] = sincos(EA.Angles(3));
    
    % Selects the matrix corresponding to the Euler rotation sequence
    if EA.Sequence == [1 2 1]
        var = [c2      0   1
                s2*s3   c3  0
                s2*c3  -s3  0];
    elseif EA.Sequence == [1 2 3]
        var = [c2*c3   s3  0
               -c2*s3   c3  0
                s2      0   1];
    elseif EA.Sequence == [1 3 1]
        var = [c2      0   1
               -s2*c3   s3  0
                s2*s3   c3  0];
    elseif EA.Sequence == [1 3 2]
        var = [c2*c3  -s3  0
               -s2      0   1
                c2*s3   c3  0];
    elseif EA.Sequence == [2 1 2]
        var = [s2*s3   c3  0
                c2      0   1
               -s2*c3   s3  0];
    elseif EA.Sequence == [2 1 3]
        var = [c2*s3   c3  0
                c2*c3  -s3  0
               -s2      0   1];
    elseif EA.Sequence == [2 3 1]
        var = [s2      0   1
                c2*c3   s3  0
               -c2*s3   c3  0];
    elseif EA.Sequence == [2 3 2]
        var = [s2*c3  -s3  0
                c2      0   1
                s2*s3   c3  0];
    elseif EA.Sequence == [3 1 2]
        var = [-c2*s3  c3  0
                s2      0   1
                c2*c3   s3  0];
    elseif EA.Sequence == [3 1 3]
        var = [s3*s2   c3  0
                s2*c3  -s3  0
                c2      0   1];
    elseif EA.Sequence == [3 2 1]
        var = [-s2     0   1
                c2*s3   c3  0
                c2*c3  -s3  0];
    elseif EA.Sequence == [3 2 3]
        var = [-s2*c3  s3  0
                s2*s3   c3  0
                c2      0   1];
    else
        error("The rotation sequence is not valid.")
    end
end

% Computes the Euler kinematic equation to transform 
% Body rotation rates to Euler Angle rates
% EA: EA type
function var = om2euler(EA)

    % Compute the sines and cosines.
    [s2, c2] = sincos(EA.Angles(2));
    [s3, c3] = sincos(EA.Angles(3));

 % Selects the matrix corresponding to the Euler rotation sequence
    if EA.Sequence == [1 2 1]
        var = [0   s3      c3
                0   s2*c3  -s2*s3
                s2 -c2*s3  -c2*c3]./s2;
    elseif EA.Sequence == [1 2 3]
        var = [c3     -s3      0
                c2*s3   c2*c3   0
               -s2*c3   s2*s3   c2]./c2;
    elseif EA.Sequence == [1 3 1]
        var = [0  -c3      s3
                0   s2*s3   s2*c3
                s2  c2*c3  -c2*s3]./s2;
    elseif EA.Sequence == [1 3 2]
        var = [c3      0   s3
               -c2*s3   0   c2*c3
                s2*c3   c2  s2*s3]./c2;
    elseif EA.Sequence == [2 1 2]
        var = [s3      0  -c3
                s2*c3   0   s2*s3
               -c2*s3   s2  c2*c3]./s2;
    elseif EA.Sequence == [2 1 3]
        var = [s3      c3      0
                c2*c3  -c2*s3   0
                s2*s3   s2*c3   c2]./c2;
    elseif EA.Sequence == [2 3 1]
        var = [0   c3     -s3
                0   c2*s3   c2*c3
                c2 -s2*c3   s2*s3]./c2;
    elseif EA.Sequence == [2 3 2]
        var = [c3      0   s3
               -s2*s3   0   s2*c3
               -c2*c3   s2 -c2*s3]./s2;
    elseif EA.Sequence == [3 1 2]
        var = [-s3     0   c3
                c2*c3   0   c2*s3
                s2*s3   c2 -s2*c3]./c2;
    elseif EA.Sequence == [3 1 3]
        var = [s3      c3      0
                s2*c3   s2*s3   0
               -c2*s3  -c2*c3   s2]./s2;
    elseif EA.Sequence == [3 2 1]
        var = [0   s3      c3
                0   c2*c3  -c2*s3
                c2  s2*s3   s2*c3]./c2;
    elseif EA.Sequence == [3 2 3]
        var = [-c3     s3      0
                s2*s3   s2*c3   0
                c2*c3  -c2*s3   s2]./s2;
    else
        error("The rotation sequence is not valid.")
    end
end

% Computes sine and cosine for improved writing and reading
% Angle: Float  (Scalar/Vector/Matrix)
function [s, c] = sincos(Angle)
    s = sin(Angle);
    c = cos(Angle);
end