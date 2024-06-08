function [p_opt, Xrank] = extract_p_opt(X)
    if any(isnan(X(:))) % load output with nan if X is nan
        p_opt = nan;
        Xrank = nan;
        return
    end
    [V,D] = eigs(X); % eigen value decomposition of X
    % sort eigenvalue/eigenvector in descending order
    Ddiag = diag(D);
    Ddiag(isnan(Ddiag)) = 0; % remove possible nan
    [d,ind] = sort(Ddiag, 'descend' );
    d = d / max(d);
    Xrank = nnz(d > 0.1);
    Dsort = D(ind,ind);
    Vsort = V(:,ind);
    % p_opt is the eigenvector corresponding to the largest eigenvalue
    v_max = Vsort(:,1);
    rho_max = Dsort(1,1);
    v_max = sqrt(rho_max)*v_max;
    t = v_max(end);
    p_opt = v_max(1:end-1); % remove the end |t|^2=1
    p_opt = p_opt / t; % correct optimal current if t has phase term
end