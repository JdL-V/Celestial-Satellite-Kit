function var = TRIAD(v1b, v2b, v1n, v2n)
    TP = types;
    bt = triadbasegen(v1b, v2b);
    nt = triadbasegen(v1n, v2n);
    var = TP.DCM(bt*nt', "B", "N");
end

function var = triadbasegen(v1, v2)
    t = zeros(3,3);
    t(:,1) = v1;
    t(:,2) = normalize(cross(v1, v2));
    t(:,3) = cross(t(:,1), t(:,2));
    var = t;
end