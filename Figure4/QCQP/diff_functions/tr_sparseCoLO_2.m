function x = tr_sparseCoLO_2(A, b, c, K, J)

    % SparseCoLO + SeDuMi
    parCoLO.SDPsolver = 'mosek';
    parCoLO.domain = 1; parCoLO.range = 2; parCoLO.EQorLMI = 1;

    [x, y, infoCoLO, cliqueDomain, cliqueRange, LOP] = sparseCoLO_Mosek(A, b, c, K, J, parCoLO);
    
%     [primalObjValue, dualObjValue, primalfeasiblity, dualfeasibility] = ...
%         evaluateCoLO(x,y,A,b,c,K,J,cliqueDomain,cliqueRange);
%     Obj = primalObjValue;

    %x = reshape(x,sqrt(length(x)),sqrt(length(x)));
    %Maybe see if have close to rank-1 when have "not PSD"
    
    x = psdCompletion(x, K, cliqueDomain);

end