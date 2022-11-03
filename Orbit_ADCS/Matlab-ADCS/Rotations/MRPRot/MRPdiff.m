function [tout, uout] = MRPdiff(f,t,u0,optiondop)
    global Dose
    Dose = @Condition;
    switch nargin
    case 4
        [tout, uout] = dosed853(f,t,u0,optiondop);
    case 3
        [tout, uout] = dosed853(f,t,u0);
    end

    function y = Condition(y)
        if dot(y(1:3), y(1:3)) > 1
            y(1:3) = -y(1:3)./dot(y(1:3),y(1:3));
        end
    end
end