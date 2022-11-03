function var = CRP2dcm(qcr)
    dcm = DCM([ 1 + qcr.x(1)^2 - qcr.x(2)^2 - qcr.x(3)^2    2*(qcr.x(1)*qcr.x(2) + qcr.x(3))            2*(qcr.x(1)*qcr.x(3) - qcr.x(2))
                2*(qcr.x(1)*qcr.x(2) - qcr.x(3))            1 - qcr.x(1)^2 + qcr.x(2)^2 - qcr.x(3)^2    2*(qcr.x(2)*qcr.x(3) + qcr.x(1)) 
                2*(qcr.x(1)*qcr.x(3) + qcr.x(2))            2*(qcr.x(2)*qcr.x(3) - qcr.x(1))            1 - qcr.x(1)^2 - qcr.x(2)^2 + qcr.x(3)^2]./(1 + dot(qcr.x,qcr.x)), ...
            qcr.outFrame, qcr.inFrame);
    var = dcm;
    % var = ((1 - dot(qcr,qcr))*Matrix(I, 3, 3) + 2*qcr*qcr' - 2*skewsym(qcr))./(1 + dot(qcr,qcr))
end