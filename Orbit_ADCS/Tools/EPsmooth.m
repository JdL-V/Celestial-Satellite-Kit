function q_array = EPsmooth(q_array, tol)

    if nargin < 2
        tol = 1;
    end

    for i = 1 : size(q_array,2) - 1
        if max(abs(q_array(:, i+1) - q_array(:, i))) > tol
            q_array(:, i+1 : end) = -q_array(:, i+1 : end);
        end
    end
end
