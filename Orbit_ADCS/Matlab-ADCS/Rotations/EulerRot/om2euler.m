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