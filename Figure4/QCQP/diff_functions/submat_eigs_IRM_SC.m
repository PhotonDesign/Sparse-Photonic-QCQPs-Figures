function Vk = submat_eigs_IRM_SC(X,bags)
    %Difference in sparseCoLO is that vectors x and xT are built into variable X
    Vk = {};
    for j = 1:length(bags)
        cur_bag = [bags{j}, length(X)]; %incorporate x-vector elements
        submat = full(X(cur_bag,cur_bag));
        [V,D] = eig(submat);
        D = diag(D);
        dim = length(D);
        Vk{j} = V(1:dim,1:(dim-1));
%         def = tw+1 - length(bags{j});
%         if def > 0
%             for k = 1:def
%                 Vk{j} = [Vk{j},zeros(size(Vk{j},1),1)];
%             end
%         end
        clear V D
    end
end