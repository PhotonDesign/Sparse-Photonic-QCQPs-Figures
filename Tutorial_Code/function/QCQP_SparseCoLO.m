function [Obj, x] = QCQP_SparseCoLO(S0, b0, c0, S, b, c, solvername)
    global timeit
    m = length(b);
    n = length(b0);
    
    % convert QCQP to a SDP
    [F0, F, t] = QCQP_2_SDP(S0, b0, c0, S, b, c);
    
    % SDP -> Conic
    n = length(b0);
    f0 = reshape(F0, (n+1)^2, 1);
    f = cell(1, m+1);
    for i = 1:m+1
        f{i} = reshape(F{i}, (n+1)^2, 1);
    end

    A = sparse(m+1, (n+1)^2); 
    for i = 1:m+1
        A(i,:) = f{i}.'; % this line needs further acceleration
    end

    b = sparse(m+1,1); b(end) = 1;
    c = f0;
    K.s = n + 1;
    J.f = m + 1;

    % SparseCoLO + SeDuMi
    parCoLO.SDPsolver = solvername;
    parCoLO.domain = 1; parCoLO.range = 2; parCoLO.EQorLMI = 1;

    timeit.init = toc;
    % [x, y, infoCoLO, cliqueDomain, cliqueRange, LOP] = sparseCoLO_Mosek(A, b, c, K, J, parCoLO);
    [x, y, infoCoLO, cliqueDomain, cliqueRange, LOP] = sparseCoLO(A, b, c, K, J, parCoLO);
    
    [primalObjValue, dualObjValue, primalfeasiblity, dualfeasibility] = ...
        evaluateCoLO(x,y,A,b,c,K,J,cliqueDomain,cliqueRange);
    Obj = primalObjValue;
    
    x = psdCompletion(x, K, cliqueDomain);
    timeit.SparseCoLO = toc - timeit.init;
end