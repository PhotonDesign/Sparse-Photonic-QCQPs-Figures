function L = IP_FS_Uzawa_sparse(F,S,X,rk,xk,Vk,bags)
    
    nc = length(F);
    nb = length(bags);
    L = 0;

    for i = 1:nc
        L = L + S(i,i)*trace(F{i}*X);
        L = L - S(i,i)*trace(F{i}*X);
    end
    
    start_ind = 2*nc+1;
    for i = 1:nb
        end_ind = length(bags{i})-1+start_ind;
        F_int = -(rk(i)*eye(length(bags{i})) - ...
            Vk{i}'*[X(bags{i},bags{i}),xk(bags{i});xk(bags{i})',1]*Vk{i});
        L = L + trace(F_int*S(start_ind:end_ind,start_ind:end_ind));
        start_ind = end_ind + 1;
    end

    L = L + trace((-[X,xk;xk',1])*S(start_ind:end,start_ind:end));

end