function [tout, uout] = EPdiff(f,t,u0,optiondop)
    global Dose
    Dose = @Condition;
    switch nargin
    case 4
        [tout, uout] = dosed853(f,t,u0,optiondop);
    case 3
        [tout, uout] = dosed853(f,t,u0);
    end

    function y = Condition(y)
        if norm(y(1:4)) ~= 1.
            y(1:4) = normalize(y(1:4), 'norm');
        end
    end
end