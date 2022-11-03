% Converts Euler angles to direction cosine matrix
% EA: EA type
function var = euler2dcm(EA)
    Mat = eye(3,3);
    for ind = 1:length(EA.Sequence)
        % Multiply the rotation matrices
        Mat = RotMat(EA.Angles(ind), EA.Sequence(ind))*Mat;
    end
    % Save the resulting matrix as a DCM type
    var = DCM(Mat, EA.outFrame, EA.inFrame);
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