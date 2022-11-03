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