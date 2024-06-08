function [Obj, x, F] = complexQCQP_SparseCoLO_MC(F0, F, t, solvername, max_or_min)
    
    if strcmp(max_or_min, 'max')
        F0 = -F0;
    end
    
    % complex -> real QCQP in dimension 2n
%     [S0, b0, c0] = c2r_QCQP(S0, b0, c0);
%     m = length(F);
%     for i = 1:m
%         [S{i}, b{i}, c{i}] = c2r_QCQP(S{i}, b{i}, c{i});
%     end
    [Obj, x, F] = QCQP_SparseCoLO2_MC(F0, F, t, solvername);
    
    if strcmp(max_or_min, 'max')
        Obj = -Obj;
    end
end