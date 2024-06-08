function [opt,rk,Xk] = IRM_Uzawa(F0,F,wk,Xk)
    F = F(1:(length(F)-1));
    dimX = length(F0);
    nc = length(F);

    [V0,D0] = eig(full(Xk));
    D0 = diag(D0);
    dim = length(D0);
    Vk = V0(1:dim,1:(dim-1));
    clear V0 D0

    S = 1/2*real(Fx_Uzawa(F,Xk,0,Vk));

    h = 1/4;

    S_last = sparse(1,1,1,2*nc+dimX-1,2*nc+dimX-1);
    %Dummy variable to start

    while norm(full(S) - full(S_last)) > 1e-7
        
       cvx_begin sdp quiet
        variable Xk(dim,dim) symmetric
        variable rk nonnegative
        minimize(trace(F0*Xk) + wk*rk + IP_FS_Uzawa(F,S,Xk,rk,Vk));
        subject to
            Xk >= 0;
            Xk(end,end) == 1;
        cvx_end
        opt = cvx_optval;

        [V0,D0] = eig(full(Xk));
        D0 = diag(D0);
        dim = length(D0);
        Vk = V0(1:dim,1:(dim-1));
        clear V0 D0

        S_last = S;
        clear S

        S = S_last + h*Fx_Uzawa(F,Xk,rk,Vk);
        [V,D] = eig(full(S));
        D(D < 0) = 0;
        D(imag(D) > 0) = 0;
        S = sparse(real(V*D*V'));
        clear V D

        h = h/2;

        %d_gap = norm(full(S) - full(S_last))
        
    end
    
end