function Fx = Fx_Uzawa(F,Xk,rk,Vk)
    dimX = size(Xk,1);
    nc = length(F);
    Fx = zeros(2*nc+dimX-1);

    for i = 1:nc
        Fx(2*i-1,2*i-1) = trace(F{i}*Xk);
        Fx(2*i,2*i) = -trace(F{i}*Xk);
    end
    
    Fx((2*nc+1):end,(2*nc+1):end) = -(rk*eye(dimX-1) - Vk'*Xk*Vk);
    
    Fx = sparse(Fx);

end