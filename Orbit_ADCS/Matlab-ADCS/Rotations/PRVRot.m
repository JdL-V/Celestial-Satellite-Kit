function var = PRVRot
    global TP 
    TP = types;

    var.PRV2dcm   = @PRV2dcm;
    var.dcm2PRV   = @dcm2PRV;
    var.PRV2om    = @PRV2om;
    var.om2PRV    = @om2PRV;
end

function var = PRV2dcm(prv)
    global TP 
    sig = 1 - cos(prv.Angle);
    dcm = TP.DCM(zeros(3,3), prv.outFrame, prv.inFrame);
    % needs revision
    for i = 1:3
        for j = 1:3
            dcm.Mat(i,j) = prv.x(i)*prv.x(j)*sig;
            if i == j
                dcm.Mat(i,j) = dcm.Mat(i,j) + cos(prv.Angle);
            else
                % dcm.Mat(i,j) = dcm.Mat(i,j) + (-1)^(i + j + (i < j))*sin(prv.Angle)*prv.x(((1 2 3)*prod.(eachcol((1 2 3).!=(i; j))))(1));
            end
        end
    end
    var = dcm;
end

function var = dcm2PRV(dcm)
    global TP 
    Angle = acos(0.5*(dcm.Mat(1,1) + dcm.Mat(2,2) + dcm.Mat(3,3) - 1));
    x = 1/(2*sin(Angle)).*[dcm.Mat(2,3)-dcm.Mat(3,2), dcm.Mat(3,1)-dcm.Mat(1,3), dcm.Mat(1,2)-dcm.Mat(2,1)];
    prv = TP.PRV(Angle, x, dcm.outFrame, dcm.inFrame);
    var = prv;
end

function var = PRV2om(prv)
    tau = skewsym(prv.x*prv.Angle);
    var = (eye(3,3) + 0.5.*tau + 1/prv.Angle^2*(1 - prv.Angle/2*cot(prv.Angle/2)).*tau^2);
end

function var = om2PRV(prv)
    tau = skewsym(prv.x.*prv.Angle);
    var = (eye(3,3) - ((1 - cos(prv.Angle))/prv.Angle^2).*tau + ((prv.Angle - sin(prv.Angle))/prv.Angle^3).*tau^2);
end