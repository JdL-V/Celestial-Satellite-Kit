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