function [Obj, x] = ...
    complexQCQP_SparseCoLO_rankPenalty(S0, b0, c0, S, b, c,...
    solvername, max_or_min,e)
    
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
    [Obj, x] = QCQP_SparseCoLO2_rankPenalty(S0, b0, c0, S, b, c, solvername,e);
    
    if strcmp(max_or_min, 'max')
        Obj = -Obj;
    end
end