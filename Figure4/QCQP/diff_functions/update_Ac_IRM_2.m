function [A,c] = update_Ac_IRM_2(wk,Xk,A,c,J)
    
    c(1) = wk;
    
    [V0,D0] = eig(full(Xk));
    D0 = diag(D0);
    dim = length(D0);
    Vk = V0(1:dim,1:(dim-1));
    
    [ncol,nrow] = size(Vk); %Relationship flips for A
    Asl = reshape(speye(nrow),nrow^2,1);
    Ass = zeros(nrow^2,ncol^2);
    for i = 1:nrow
        for j = 1:nrow
            row = nrow*(i-1)+j;
            Ass(row,:) = -reshape(Vk(:,i)*Vk(:,j)',1,ncol^2);
        end
    end
    A((J.f+1):end,:) = [Asl,Ass];
    A = sparse(A);
    c = sparse(c);
end