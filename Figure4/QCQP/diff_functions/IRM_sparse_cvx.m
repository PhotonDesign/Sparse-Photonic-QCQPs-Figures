function [opt,rk,Xk,xk] = IRM_sparse_cvx(F0,F,wk,Vk,bags)
    dim = length(F0);
    cvx_begin sdp quiet
        variable Xk(dim,dim) symmetric
        variable rk(length(bags)) nonnegative
        variable xk(dim)
        minimize(trace(F0*Xk) + wk*norm(rk,1));
        subject to
            for i = 1:(length(F)-1) %Put in last constraint manually
                trace(F{i}*Xk) == 0;
            end
            %Xk >= 0;
            Xk(end,end) == 1;
            [Xk,xk;xk',1] >= 0; %Takes care of 2 commented out sets of constraints
            for j = 1:length(bags)
                %[Xk(bags{j},bags{j}),xk(bags{j});xk(bags{j})',1] >= 0;
                rk(j)*eye(length(bags{j})) -...
                    Vk{j}'*[Xk(bags{j},bags{j}),xk(bags{j});xk(bags{j})',1]*Vk{j} >= 0;
            end
        cvx_end
    opt = cvx_optval;
end
