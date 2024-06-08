function [Obj, x] = complexQCQP_SparseCoLO(S0, b0, c0, S, b, c, solvername, max_or_min)
% function [Obj] = complexQCQP_SparseCoLO(S0, b0, c0, S, b, c, solvername, max_or_min)
    
    if strcmp(max_or_min, 'max')
        S0 = -S0;
        b0 = -b0;
        c0 = -c0;
    end
    
    % complex -> real QCQP in dimension 2n
    [S0, b0, c0] = c2r_QCQP(S0, b0, c0);
    m = length(S);
    for i = 1:m
        [S{i}, b{i}, c{i}] = c2r_QCQP(S{i}, b{i}, c{i});
    end
    
    
    [Obj, x] = QCQP_SparseCoLO(S0, b0, c0, S, b, c, solvername);
%     [Obj] = QCQP_SparseCoLO(S0, b0, c0, S, b, c, solvername);
    


    if strcmp(max_or_min, 'max')
        Obj = -Obj;
    end
end