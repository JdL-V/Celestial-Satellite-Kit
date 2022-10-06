function var = getCM(M, R)
    % Works for position, speed and acceleration
    CM = zeros(3,1);
    for j = 1:length(M)
        CM = CM + M(j).*R(:,j);
    end
    var = CM/sum(M);
end

function var = getForce(M, ddRc)
    var = sum(M).*ddRc;
end

function var = getKineticE(M, dR)
    dRc = getCM(M, dR)
    T1 = 0.5*sum(M)*dot(dRc, dRc)

    T2 = 0
    T = 0
    for j = 1:length(M)
        % T = T + 0.5*M(j)*dot(dR(j), dR(j))
        T2 = T2 + 0.5*M(j)*dot(dR(j) - dRc, dR(j) - dRc)
    end
    T = T1 + T2
    var = [T, T1, T2]
end

function var = getLinearMom(M, dR)
    dRc = getCM(M, dR)
    var = sum(M).*dRc
end

function var = getAngularMom(M, R, dR, Rp, dRp)
    H = zeros(3)
    for j = 1:length(M)
        H = H + M(j).*cross(R(j) - Rp, dR(j) - dRp)
    end
    var = H
end

function var = getTorque(M, R, Rp, F)
    Lp = zeros(3,1);
    for j = 1:length(M)
        Lp = Lp + cross(R(j) - Rp, F);
    end
    var = Lp;
end

% function var = getInertialtDev(M, R, Rp, F)
%     Lp = zeros(3)
%     for j = 1:length(M)
%         Lp = Lp + cross(R(j) - Rp, F)
%     end
%     var = Lp
% end