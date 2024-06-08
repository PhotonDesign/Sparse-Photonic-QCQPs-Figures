function [fobj, popt, Xrank] = QCQP_solver(S0,b0,c0,S,b,c)
    % Solve a nonconvex QCQP via semidefinite relaxation
    % Maximize    x'*S0*x + Re(b0'*x) + c0
    % subject to  x'*Si*x + Re(bi'*x) + ci = 0, i = 1, ..., N
    % Si (i = 0, 1, ..., N) are Hermittian matrices

    H0 = [S0, b0/2; b0'/2 c0];
    NS = length(S);
    for i = 1:NS
        H{i} = [S{i}, b{i}/2; b{i}'/2 c{i}];
    end
    
    N = length(b0);
    
    cvx_begin sdp % quiet 
        variable X(N+1,N+1) hermitian
        maximize( real(H0(:)'*X(:)) );
        subject to
            for i = 1:NS
                H{i}(:)'*X(:) == 0;
            end
            X(N+1,N+1) == 1;
            X >= 0;
    cvx_end
    fobj = cvx_optval;
    
    [popt, Xrank] = extract_p_opt(X);
end

