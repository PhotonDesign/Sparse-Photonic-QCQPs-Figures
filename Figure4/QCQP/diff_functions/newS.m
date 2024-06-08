function S = newS(S0,F,nc,bags,Xk,rk,xk,Vk,h)
    S = zeros(size(S0));
    for i = 1:nc
        S(2*i-1,2*i-1) = S0(2*i-1,2*i-1) + h*(trace(F{i}*Xk));
        S(2*i,2*i) = S0(2*i,2*i) + h*(- trace(F{i}*Xk));
    end
    start_ind = 2*nc+1;
    for i = 1:length(bags)
        end_ind = length(bags{i})-1+start_ind; 
        S(start_ind:end_ind,start_ind:end_ind) =...
 S0(start_ind:end_ind,start_ind:end_ind) + h*(-(rk(i)*eye(length(bags{i})) - ...
            Vk{i}'*[Xk(bags{i},bags{i}),xk(bags{i});xk(bags{i})',1]*Vk{i}));
        start_ind = end_ind + 1;
    end
    S(start_ind:end,start_ind:end) = ...
S0(start_ind:end,start_ind:end) + h*(-[Xk,xk;xk',1]);

    [V,D] = eig(S);
    D(D < 0) = 0;
    D(imag(D) > 0) = 0;
    S = sparse(V*D*V');

end