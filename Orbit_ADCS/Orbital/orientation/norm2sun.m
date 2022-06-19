function Panel = norm2sun()
    npanels = 3;
    Panel = zeros(size(th,2),npanels);
    for k = 1:npoints
        irr = a_e^2*irr0/norm(R(k,:))^2;
        Panel(k,:) = ones(1,npanels)*irr*visibility(k);
    end