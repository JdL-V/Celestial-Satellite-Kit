% This modified function computes atan with the result depending on the quadrant
% y: Float
% x: Float
function var = atanMk2(y, x) 
    var = wrapTo2Pi(atan((y)/(x)) + pi*(sign(x) == sign(y))*(sign(y) == -1) ...
                                  - pi*(sign(x) ~= sign(y))*(sign(y) == -1));
end