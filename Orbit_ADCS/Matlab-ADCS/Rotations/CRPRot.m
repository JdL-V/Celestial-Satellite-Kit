function var = CRPRot
    global TP 
    TP = types;

    var.CRP2dcm   = @CRP2dcm;
    var.dcm2CRP   = @dcm2CRP;
    var.CRP2PRV   = @CRP2PRV;
    var.PRV2CRP   = @PRV2CRP;
    var.SumCRP    = @SumCRP;
    var.CRP2om    = @CRP2om;
    var.om2CRP    = @om2CRP;
end

function var = CRP2dcm(qcr)
    global TP 
    dcm = TP.DCM([ 1 + qcr.x(1)^2 - qcr.x(2)^2 - qcr.x(3)^2    2*(qcr.x(1)*qcr.x(2) + qcr.x(3))            2*(qcr.x(1)*qcr.x(3) - qcr.x(2))
                2*(qcr.x(1)*qcr.x(2) - qcr.x(3))            1 - qcr.x(1)^2 + qcr.x(2)^2 - qcr.x(3)^2    2*(qcr.x(2)*qcr.x(3) + qcr.x(1)) 
                2*(qcr.x(1)*qcr.x(3) + qcr.x(2))            2*(qcr.x(2)*qcr.x(3) - qcr.x(1))            1 - qcr.x(1)^2 - qcr.x(2)^2 + qcr.x(3)^2]./(1 + dot(qcr.x,qcr.x)), ...
            qcr.outFrame, qcr.inFrame);
    var = dcm;
    % var = ((1 - dot(qcr,qcr))*Matrix(I, 3, 3) + 2*qcr*qcr' - 2*skewsym(qcr))./(1 + dot(qcr,qcr))
end

function var = dcm2CRP(dcm)
    EP = EPRot;
    % ζ = sqrt(1 + dcm[1,1] + dcm[2,2] + dcm[3,3])
    % var = [dcm[2,3] - dcm[3,2]
    %         dcm[3,1] - dcm[1,3] 
    %         dcm[1,2] - dcm[2,1]]./ζ^2
    var = EP2CRP(EP.dcm2EP(dcm));
end

function var = CRP2EP(qcr)
    global TP 
    q = TP.quaternion([1, qcr.x(1:3)]./sqrt(1 + dot(qcr,qcr)), ...
                    qcr.outFrame, qcr.inFrame);
    var = q;
end

function var = EP2CRP(q)
    global TP 
    qcr = TP.CRP(q.x(2:4)./q.x(1), q.outFrame, q.inFrame);
    var = qcr;
end

function var = CRP2PRV(qcr)
    global TP 
    T = atan(qcr.x);
    e = normalize(T, 'norm');
    prv = TP.PRV(2*T(1)/e(1), e, qcr.outFrame, qcr.inFrame);
    var = prv;
end

function var = PRV2CRP(prv)
    global TP 
    qcr = TP.CRP(tan(prv.Angle/2).*prv.x, prv.outFrame, prv.inFrame);
    var = qcr;
end

function var = SumCRP(qcr1, qcr2)
    global TP 
    qcr = TP.CRP((qcr2.x + qcr1.x - cross(qcr2.x,qcr1.x))./(1 - dot(qcr2.x,qcr1.x)), ...
                qcr2.outFrame, qcr1.inFrame);
    var = qcr;
end

function var = CRP2om(qcr)
    var = 2/(1 + dot(qcr, qcr)).*(eye(3,3) - skewsym(qcr));
end

function var = om2CRP(qcr)
    var = [1 + qcr(1)^2            qcr(1)*qcr(2) - qcr(3)  qcr(1)*qcr(3) + qcr(2)
            qcr(2)*qcr(1) + qcr(3)  1 + qcr(2)^2            qcr(2)*qcr(3) - qcr(1)
            qcr(3)*qcr(1) - qcr(2)  qcr(3)*qcr(2) + qcr(1)  1 + qcr(3)^2          ].*0.5;
    % var = 0.5.*(Matrix(I, 3, 3) + skewsym(qcr) + qcr*qcr')
end