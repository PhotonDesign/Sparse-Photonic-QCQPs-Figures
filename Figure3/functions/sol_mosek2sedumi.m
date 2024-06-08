function [x, y, obj] = sol_mosek2sedumi(res, Ks)
    
    % extract optimal y and obj
    y = res.sol.itr.y;
    obj = res.sol.itr.pobjval;
    
    % extract optimal x in SparseCoLO format
    barx = res.sol.itr.barx;
    Nk = length(Ks);
    idx1 = 1;
    x = [];
    for i = 1:Nk
        k = Ks(i);
        d = (k+1)*k/2;
        idx2 = idx1 + d - 1;
        barx_i = barx(idx1:idx2);
        idx1 = idx1 + d;

        X = zeros(k);
        X(tril(true(k))) = barx_i;
        X = X + tril(X,-1)';
        x = [x; X(:)];
    end
end