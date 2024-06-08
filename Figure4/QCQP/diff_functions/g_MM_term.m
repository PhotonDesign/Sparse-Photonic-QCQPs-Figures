function g = g_MM_term(Wk,eps,eta)
    [P,D] = eig(full(Wk));
    %D(D < 0) = 0; %Get rid of any trivially small non-PSD elements
    g = eta*(1/eps)*P*diag(exp(-diag(D./eps)))*P';
end