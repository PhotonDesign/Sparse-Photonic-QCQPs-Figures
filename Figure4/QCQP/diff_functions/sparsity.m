function [i,j,x] = sparsity(X)
    %Get sparsity of a cell of matrix constraints
    if iscell(X)
        x = sparse(zeros(size(X{1})));
        for k = 1:length(X)
            x = x + X{k};
        end
        [i,j] = find(x);
    else
        [i,j] = find(X);
        x = X;
    end
end