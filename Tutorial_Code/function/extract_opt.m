function [es_opt, d, X_opt] = extract_opt(x_opt, Nx)
    X_opt = reshape(x_opt, 2*Nx+1, 2*Nx+1);
    [V, D] = eigs(X_opt);
    d = diag(D);
    x1 = V(:,1);
    x1 = x1 * sqrt(d(1));
    x1 = x1 / x1(end);
    x1_complex = x1(1:Nx) + 1j*x1(Nx+1:2*Nx);
    es_opt = x1_complex;
end