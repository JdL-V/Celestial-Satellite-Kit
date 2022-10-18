function var = MRPRot
    global TP 
    TP = types;

    var.MRP2dcm   = @MRP2dcm;
    var.dcm2MRP   = @dcm2MRP;
    var.MRP2PRV   = @MRP2PRV;
    var.PRV2MRP   = @PRV2MRP;
    var.MRP2CRP   = @MRP2CRP;
    var.CRP2MRP   = @CRP2MRP;
    var.SumMRP    = @SumMRP;
    var.MRP2om    = @MRP2om;
    var.om2MRP    = @om2MRP;
    var.MRPdiff   = @MRPdiff;
end

function var = MRP2dcm(qmr)
    global TP 
    qmrMat = skewsym(qmr.x);
    qmrSq = dot(qmr.x,qmr.x);
    dcm = TP.DCM(eye(3, 3) + (8*qmrMat^2 - 4*(1 - qmrSq)*qmrMat)/(1 + qmrSq)^2, ...
            qmr.outFrame, qmr.inFrame);
    var = dcm;
end

function var = dcm2MRP(dcm)
    EP = EPRot;
    % ζ = sqrt(1 + dcm[1,1] + dcm[2,2] + dcm[3,3])
    % var = [dcm[2,3] - dcm[3,2]
    %         dcm[3,1] - dcm[1,3] 
    %         dcm[1,2] - dcm[2,1]]./(ζ*(ζ + 2))
    var = EP2MRP(EP.dcm2EP(dcm));
end

function var = MRP2EP(qmr)
    global TP 
    sig = dot(qmr.x,qmr.x);
    q = TP.quaternion([(1 - sig)/(1 + sig) 2*qmr(1:3)'/(1 + sig)], ...
                    qmr.outFrame, qmr.inFrame);
    var = q;
end

function var = EP2MRP(q)
    global TP 
    qmr = TP.MRP(q.x(2:4)./(1 + q.x(1)), ...
                q.outFrame, q.inFrame);
    var = qmr;
end

function var = MRP2CRP(qmr)
    global TP 
    qcr = TP.CRP(qmr.x./(1 - dot(qmr.x,qmr.x)), ...
                qmr.outFrame, qmr.inFrame);
    var = qcr;
end

function var = CRP2MRP(qcr)
    global TP 
    qmr = TP.MRP(qcr.x./(1 + sqrt(1 + dot(qcr.x,qcr.x))), ...
                qcr.outFrame, qcr.inFrame);
    var = qmr;
end

function var = MRP2PRV(qmr)
    global TP 
    T = atan(qmr.x);
    e = normalize(T, 'norm');
    prv = TP.PRV(4*T(1)/e(1), e, qmr.outFrame, qmr.inFrame);
    var = prv;
end

function var = PRV2MRP(prv)
    global TP 
    qmr = TP.MRP(tan(prv.Angle/4).*prv.x, ...
                prv.outFrame, prv.inFrame);
    var = qmr;
end

function var = SumMRP(qmr1, qmr2)
    global TP 
    qmr1Sq = dot(qmr1.x,qmr1.x);
    qmr2Sq = dot(qmr2.x,qmr2.x);
    qmr = TP.MRP(((1 - qmr1Sq)*qmr2.x + (1 - qmr2Sq)*qmr1.x - 2*cross(qmr2.x,qmr1.x)) ...
            ./(1 + qmr1Sq*qmr2Sq - 2*dot(qmr1.x,qmr2.x)), ...
            qmr2.outFrame, qmr1.inFrame);
    var = qmr;
end

function var = MRP2om(qmr)
    var = (4/(1 + dot(qmr,qmr))^2)*(om2MRP(qmr)');
end

function var = om2MRP(qmr)
    qmrSq = dot(qmr,qmr);
    var = [1 - qmrSq + 2*qmr(1)^2      2*(qmr(1)*qmr(2) - qmr(3))  2*(qmr(1)*qmr(3) + qmr(2))
            2*(qmr(2)*qmr(1) + qmr(3))  1 - qmrSq + 2*qmr(2)^2      2*(qmr(2)*qmr(3) - qmr(1))
            2*(qmr(3)*qmr(1) - qmr(2))  2*(qmr(3)*qmr(2) + qmr(1))  1 - qmrSq + 2*qmr(3)^2    ]./4;
    % var = ((1 - qmrSq).*Matrix(I,3,3) + 2*skewsym(qmr) + 2*qmr*qmr')./4
end

% function var = MRPdiff(nn::Int64, tmax, u0, f::Function var =, sch::Function var =)
%     h  = tmax/(nn-1)
%     th::LinRange{Float64, Int64} = LinRange(0, tmax, nn)
%     u  = zeros(nn, length(u0))
%     u[1,:] = u0

%     for n::Int64 = 1:(nn-1)
%         u[n+1,:] = sch(h, u[n,:], th[n], f)
%         % Avoid singularity (1:3 in order to use the equations with omega rates)
%         if dot(u[n+1,1:3], u[n+1,1:3]) > 1
%             u[n+1,1:3] = -u[n+1,1:3]./dot(u[n+1,1:3],u[n+1,1:3])
%         end
%     end
%     var = th, u
% end

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