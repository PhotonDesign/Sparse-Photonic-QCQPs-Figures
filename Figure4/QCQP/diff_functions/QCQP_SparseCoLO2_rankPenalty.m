function [Obj, x] = ...
    QCQP_SparseCoLO2_rankPenalty(S0, b0, c0, S, b, c, solvername,e)
    
    m = length(b);
    n = length(b0);
    
    % convert QCQP to a SDP
    [F0, F, t] = QCQP_2_SDP(S0, b0, c0, S, b, c);
    
    %Add penalty term
    F0 = F0 - e*speye(size(F0));
    
    n = length(b0);
    f0 = reshape(F0, (n+1)^2, 1);
    f = cell(1, m+1);
    for i = 1:m+1
        f{i} = reshape(F{i}, (n+1)^2, 1);
    end
    
    row_idx_A = [];
    col_idx_A = [];
    val_A = [];
    for i = 1:m+1
        [row_idx,col_idx,val] = find(f{i}.');
        row_idx = row_idx + (i-1);
        row_idx_A = [row_idx_A, row_idx];
        col_idx_A = [col_idx_A, col_idx];
        val_A = [val_A, val];
    end
    
    A = sparse(row_idx_A, col_idx_A, val_A, m+1, (n+1)^2);

    b = sparse(m+1,1); b(end) = 1;
    c = f0;
    K.s = n + 1;
    J.f = m + 1;

    % SparseCoLO + SeDuMi
    parCoLO.SDPsolver = solvername;
    parCoLO.domain = 1; parCoLO.range = 2; parCoLO.EQorLMI = 1;

    [x, y, infoCoLO, cliqueDomain, cliqueRange, LOP] = sparseCoLO_Mosek(A, b, c, K, J, parCoLO);
    
    [primalObjValue, dualObjValue, primalfeasiblity, dualfeasibility] = ...
        evaluateCoLO(x,y,A,b,c,K,J,cliqueDomain,cliqueRange);
    Obj = primalObjValue;
    
    x = psdCompletion(x, K, cliqueDomain);
    
end
