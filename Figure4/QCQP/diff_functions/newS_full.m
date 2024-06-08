function S = newS_full(F,S,Xk,rk,Vk)
    for i = 1:nc
        S(i,i) = S(i,i)*trace(F{i}*Xk);
        S(i,i) = - S(i,i)*trace(F{i}*Xk);
    end
end