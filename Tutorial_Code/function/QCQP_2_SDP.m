function [F0, F, t] = QCQP_2_SDP(S0, b0, c0, S, b, c)
    % convert a QCQP of the following form, x in R^n or C^n
    %  Min  x'*S0*x + real(b0*x) + c0
    %  s.t. x'*S{i}*x + real(b{i}*x) + c{i} = 0,  i = 1, ..., m 
    % to a SDP of the following form
    %  Min  Tr(F0*X)
    %  s.t. Tr(F{i}*X) = t{i}, i = 1 = 1, ..., m

    m = length(S);
    n = length(b0);
    
    F0 = [S0, b0/2; b0'/2, c0]; 
    F = cell(1, m+1);
    t = cell(1, m+1);
    for i = 1:m
        F{i} = [S{i}, b{i}/2; b{i}'/2, c{i}];
        t{i} = 0;
    end
    F{m+1} = sparse(n+1, n+1, 1, n+1, n+1);
    t{m+1} = 1;
end