function [Obj, x, F] = QCQP_SparseCoLO2_MC(F0, F, t, solvername)
    m = length(F); %F is a cell array
    n = length(F0); %F0 is a square matrix
    
    % convert QCQP to a SDP
    %[F0, F, t] = QCQP_2_SDP(S0, b0, c0, S, b, c);
    
    f0 = reshape(F0, n^2, 1);
    f = cell(1, m);
    for i = 1:m
        f{i} = reshape(F{i}, n^2, 1);
    end
    
    row_idx_A = [];
    col_idx_A = [];
    val_A = [];
    for i = 1:m
        [row_idx,col_idx,val] = find(f{i}.');
        row_idx = row_idx + (i-1);
        row_idx_A = [row_idx_A, row_idx];
        col_idx_A = [col_idx_A, col_idx];
        val_A = [val_A, val];
    end
    
    A = sparse(row_idx_A, col_idx_A, val_A, m, n^2);

    b = sparse(cell2mat(t));
    c = f0;
    K.s = n;
    J.f = m;

    % SparseCoLO + SeDuMi
    parCoLO.SDPsolver = solvername;
    parCoLO.domain = 1; parCoLO.range = 2; parCoLO.EQorLMI = 1;

    [x, y, infoCoLO, cliqueDomain, cliqueRange, LOP] = sparseCoLO_Mosek(A, b, c, K, J, parCoLO);
    
    [primalObjValue, dualObjValue, primalfeasiblity, dualfeasibility] = ...
        evaluateCoLO(x,y,A,b,c,K,J,cliqueDomain,cliqueRange);
    Obj = primalObjValue;
    
    x = psdCompletion(x, K, cliqueDomain);
    
end
