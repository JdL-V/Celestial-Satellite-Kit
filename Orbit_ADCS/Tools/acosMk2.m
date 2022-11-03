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