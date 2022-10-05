function var = EPRot
    global TP 
    TP = types;

    var.EP2dcm   = @EP2dcm;
    var.dcm2EP   = @dcm2EP;
    var.EP2PRV   = @EP2PRV;
    var.PRV2EP   = @PRV2EP;
    var.SumEP    = @SumEP;
    var.EP2om    = @EP2om;
    var.om2EP    = @om2EP;
    var.EPdiff   = @EPdiff;
end

function var = EP2dcm(q)
    global TP 
    dcm = TP.DCM([q.x(1)^2 + q.x(2)^2 - q.x(3)^2 - q.x(4)^2   2*(q.x(2)*q.x(3) + q.x(1)*q.x(4))           2*(q.x(2)*q.x(4) - q.x(1)*q.x(3))
            2*(q.x(2)*q.x(3) - q.x(1)*q.x(4))           q.x(1)^2 - q.x(2)^2 + q.x(3)^2 - q.x(4)^2   2*(q.x(3)*q.x(4) + q.x(1)*q.x(2)) 
            2*(q.x(2)*q.x(4) + q.x(1)*q.x(3))           2*(q.x(3)*q.x(4) - q.x(1)*q.x(2))           q.x(1)^2 - q.x(2)^2 - q.x(3)^2 + q.x(4)^2],...
            q.outFrame, q.inFrame);
    var = dcm;
end

% function var = sheppard(dcm)
%     trace = dcm.Mat(1,1) + dcm.Mat(2,2) + dcm.Mat(3,3)
%     q = quaternion( sqrt.([ 1 + trace
%                             1 - trace + 2*dcm.Mat(1,1)
%                             1 - trace + 2*dcm.Mat(2,2)
%                             1 - trace + 2*dcm.Mat(3,3)]./4),
%                     dcm.outFrame,
%                     dcm.inFrame)

%     index = (1 2 3 4)[findfirst(q.x .== maximum(q.x))]
%     if index == 1
%         for j = 1:3
%             i = (1:3)(1:end .~= j]
%             q.x[j+1)  = (dcm.Mat[i(1),i(2)] - dcm.Mat[i(2),i(1)])/(4*q.x[index])*(-1).^(j + 1)
%         end
%     else
%         i = (1:3)(1:end .~= index-1)
%         q.x(1) = (dcm.Mat[i(1),i(2)] - dcm.Mat[i(2),i(1)])/(4*q.x[index])*(-1).^index
%         for j = (1:3)(1:end .~= index-1)
%             q.x[j+1) = (dcm.Mat[index-1,j] + dcm.Mat[j,index-1))/(4*q.x[index])
%         end
%         q.x = (q.x).*sign(q.x(1))
%     end
%     var = q
% end

function var = dcm2EP(dcm, eta)
    global TP
    trace = dcm.Mat(1,1) + dcm.Mat(2,2) + dcm.Mat(3,3);
    check = [trace,
             2*dcm.Mat(1,1) - trace,
             2*dcm.Mat(2,2) - trace,
             2*dcm.Mat(3,3) - trace];

    q = TP.quaternion(zeros(4), dcm.outFrame, dcm.inFrame);
    for i = 1:4
        k0 = 1:3;
        k = k0(1:end ~= (i-1));
        sgn = sign((dcm.Mat(k(1),k(2)) - dcm.Mat(k(2),k(1)))*(-1).^i)*(i ~= 1) + (i == 1);

        if check(i) > eta
            q.x(i) = 0.5*sqrt(1 + check(i))*sgn;
        else
            q.x(i) = 0.5*sqrt(( (dcm.Mat(3,2) + dcm.Mat(2,3)*(-1)^(i == 2 || i == 1))^2 ... 
                              + (dcm.Mat(1,2) + dcm.Mat(2,1)*(-1)^(i == 4 || i == 1))^2 ... 
                              + (dcm.Mat(3,1) + dcm.Mat(1,3)*(-1)^(i == 3 || i == 1))^2) / (3 - check(i)))*sgn;
        end
    end
    var = q;
end

function var = EP2PRV(q)
    global TP
    angle = 2*acos(q.x(1));
    prv = TP.PRV(angle, q.x(2:4)./sin(angle/2),...
                q.outFrame, q.inFrame);
    var = prv;
end

function var = PRV2EP(prv)
    global TP 
    q = TP.quaternion([cos(prv.Angle/2), prv.x(1)*sin(prv.Angle/2), prv.x(2)*sin(prv.Angle/2), prv.x(3)*sin(prv.Angle/2)],...
                        prv.outFrame, prv.inFrame);
    var = q;
end

function var = SumEP(q1, q2)
    global TP 
    q = TP.quaternion([q2.x(1) -q2.x(2) -q2.x(3) -q2.x(4);
                    q2.x(2)  q2.x(1)  q2.x(4) -q2.x(3);
                    q2.x(3) -q2.x(4)  q2.x(1)  q2.x(2);
                    q2.x(4)  q2.x(3) -q2.x(2)  q2.x(1)]*(q1.x)',...
                    q2.outFrame, q1.inFrame);
    var = q;
end

function var = EP2om(q)
    var = inv([q(1) -q(2) -q(3) -q(4)
                q(2)  q(1) -q(4)  q(3)
                q(3)  q(4)  q(1) -q(2)
                q(4) -q(3)  q(2)  q(1)]./2);
end

function var = om2EP(q)
    var = [-q(2) -q(3) -q(4)             % [q(1) -q(2) -q(3) -q(4)
             q(1) -q(4)  q(3)             %  q(2)  q(1) -q(4)  q(3)
             q(4)  q(1) -q(2)             %  q(3)  q(4)  q(1) -q(2)
            -q(3)  q(2)  q(1)]./2;        %  q(4) -q(3)  q(2)  q(1)]
end

function [tout, uout] = EPdiff(f,t,u0,optiondop)
    global Dose
    Dose = @Condition;
    switch nargin
    case 4
        [tout, uout] = dop853(f,t,u0,optiondop);
    case 3
        [tout, uout] = dop853(f,t,u0);
    end

    function y = Condition(y)
        if norm(y(1:4)) ~= 1.
            y(1:4) = normalize(y(1:4), 'norm');
        end
    end
end