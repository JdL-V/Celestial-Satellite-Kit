function [Vals,Vecs] = eigenSort(Vals, Vecs)
    Vals = diag(Vals);
    if ~issorted(Vals)
        [Vecs, Vals] = eig(A);
        [Vals, Inds] = sort(Vals);
        Vecs = Vecs(:, Inds);
    end
end