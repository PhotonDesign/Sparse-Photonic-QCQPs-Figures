function [es_opt, d, X_opt] = extract_opt_sparsity(x_opt, Nx, dim2)
    X_opt = reshape(x_opt, dim2, dim2);
    X_opt = X_opt(1:(2*Nx+1),1:(2*Nx+1)); %Get proper submatrix
    [V, D] = eigs(X_opt);
    d = diag(D);
    x1 = V(:,1);
    x1 = x1 * sqrt(d(1));
    x1 = x1 / x1(end); %%FIx this tomorrow!
    x1_complex = x1(1:Nx) + 1j*x1(Nx+1:2*Nx);
    es_opt = x1_complex;
end