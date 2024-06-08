function Vk = submat_eigs_IRM(X,x,bags)
    Vk = {};
    for j = 1:length(bags)
        submat = X(bags{j},bags{j});
        subvec = x(bags{j});
        ex_submat = [submat,subvec;subvec',1];
        [V,D] = eig(full(ex_submat));
        D = diag(D);
        dim = length(D);
        Vk{j} = V(1:dim,1:(dim-1));
        clear V D
    end
end