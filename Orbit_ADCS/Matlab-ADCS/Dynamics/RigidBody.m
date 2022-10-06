function var = getInertiaTensor(M, R)
    Io = zeros(3,3);
    for j = 1:length(M)
        Io = Io - M(j).*cross2mat(R(:,j))^2;
    end
    var = I;
end

function var = translateInertiaFrame(Ig, M, Rp)
    var = Ig + sum(M)*cross2mat(Rp)*cross2mat(Rp)';
end

function var = rotateInertiaFrame(Io, dcm)
    var = dcm*Io*dcm';
end

function var = getPrincipalInertia(Io)
    [dcm, Lam] = eig(Io)';
    [Lam, dcm] = eigenSort(Lam, dcm);
    % Right handed coordinate frame check
    dcm(3,:) = cross(dcm(1,:), dcm(2,:));
    var = diagm(Lam), dcm;
end

function var = getRotKineticE(Io, om)
    om = om(:);
    var = 0.5*om'*Io*om;
end