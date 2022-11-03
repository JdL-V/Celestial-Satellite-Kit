function [var1, var2] = getPrincipalInertia(Io)
    [dcm, Lam] = eig(Io);
    [Lam, dcm] = eigenSort(Lam, dcm, 'descend');
    dcm = dcm';
    % Right handed coordinate frame check
    dcm(3,:) = cross(dcm(1,:), dcm(2,:));
    var1 = diag(Lam);
    var2 = dcm;
end