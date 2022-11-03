function var = getInertiaTensor(M, R)
    Io = zeros(3,3);
    for j = 1:length(M)
        Io = Io - M(j).*skewsym(R(:,j))^2;
    end
    var = I;
end