function [Vals,Vecs] = eigenSort(Vals, Vecs, direction)
    if nargin < 3
        direction = 'ascend';
    end
    Vals = diag(Vals);
    [Vals, Inds] = sort(Vals, direction);
    Vecs = Vecs(:, Inds);
end