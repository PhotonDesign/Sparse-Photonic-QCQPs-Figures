function [opt,rk,Xk] = IRM(F0,F,wk,Vk)
    dim = length(F0);
    cvx_begin sdp quiet
        variable Xk(dim,dim) symmetric
        variable rk nonnegative
        minimize(trace(F0*Xk) + wk*rk);
        subject to
            for i = 1:(length(F)-1) %Put in last constraint manually
                trace(F{i}*Xk) == 0;
            end
            Xk >= 0;
            Xk(end,end) == 1;
            rk*eye(dim-1) - Vk'*Xk*Vk >= 0;
        cvx_end
    opt = cvx_optval;
end
