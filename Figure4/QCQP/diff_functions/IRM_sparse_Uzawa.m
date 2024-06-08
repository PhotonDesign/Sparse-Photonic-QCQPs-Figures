function [cvx_optval,rk,Xk,xk,S] = IRM_sparse_Uzawa(F0,F,wk,X0,bags)
    
    F = F(1:(length(F)-1));
    dimX = length(F0);
    nb = length(bags);
    nc = length(F);
    
    sb = 0;
    for i = 1:nb
        sb = sb + length(bags{i});
    end

    [V0,D0] = eigs(X0);
    x0 = V0(:,1)*sqrt(D0(1,1));
    Vk = submat_eigs_IRM(X0,x0,bags);
    rk = 100*ones(1,nb);
    
    h = 1;
    S = sparse(2*nc+sb+dimX+1,2*nc+sb+dimX+1);
    S = newS(S,F,nc,bags,X0,rk,x0,Vk,h); %Start with rk = 0

    S_last = sparse(1,1,1,2*nc+sb+dimX+1,2*nc+sb+dimX+1);
    %Dummy variable to start

  % while norm(full(S) - full(S_last)) > 1e-7 %Stopping condition: S = new S
        
        cvx_begin sdp quiet
            variable Xk(dimX,dimX) symmetric
            variable rk(nb) nonnegative
            variable xk(dimX)
            minimize(trace(F0*Xk) + wk*sum(rk) + IP_FS_Uzawa_sparse(F,S,Xk,rk,xk,Vk,bags))
            subject to
                Xk(end,end) == 1;
                Xk >= 0;
        cvx_end

        S_last = S;
        clear S

        S = newS(S_last,F,nc,bags,Xk,rk,xk,Vk,h);
        
        h = h/2;
            
        norm(full(S) - full(S_last))

   % end
    
end
