function IP = IP_FS_Uzawa(F,S,Xk,rk,Vk)
    dimX = size(Xk,1);
    nc = length(F);

    IP = 0;

    for i = 1:nc
        IP = IP + S(2*i-1,2*i-1)*trace(F{i}*Xk);
        IP = IP + S(2*i,2*i)*(-trace(F{i}*Xk));
    end

    IP = IP + trace((-(rk*eye(dimX-1) - Vk'*Xk*Vk))*S((2*nc+1):end,(2*nc+1):end));

end