function var = getCM(M, R)
    % Works for position, speed and acceleration
    CM = zeros(3,1);
    for j = 1:length(M)
        CM = CM + M(j).*R(:,j);
    end
    var = CM/sum(M);
end