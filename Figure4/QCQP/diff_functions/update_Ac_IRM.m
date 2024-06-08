function [A,c] = update_Ac_IRM(wk,Xk,bags,A,c,J,nrow_As)
    nb = length(bags);
    c(1:nb) = ones(nb,1)*wk;
    Vk = submat_eigs_IRM_SC(Xk,bags);
    Asl = sparse(nrow_As,nb);
    Ass = sparse(nrow_As,length(Xk)^2);
    row_start = 1;
    for k = 1:nb
        cur_bag = [bags{k}, length(Xk)];
        bag_size = length(cur_bag);
        noVs = bag_size - 1;
        row_end = row_start + noVs^2 -1;
        Asl(row_start:row_end,k) = reshape(speye(noVs),noVs^2,1);
        Ass_int = zeros(noVs^2,bag_size^2);
        spinds_int = zeros(bag_size^2,length(Xk)^2);
        for i = 1:noVs
            for j = 1:noVs
                row = noVs*(i-1)+j;
                Ass_int(row,:) = -reshape(Vk{k}(:,i)*Vk{k}(:,j)',bag_size^2,1);
            end
        end
        for i = 1:bag_size
            for j = 1:bag_size
                row = bag_size*(i-1)+j;
                ind = sub2ind(size(Xk),cur_bag(j),cur_bag(i));
                spinds_int(row,:) = sparse(1,ind,1,1,length(Xk)^2);
            end
        end
        Ass(row_start:row_end,:) = Ass_int*spinds_int;
        clear Ass_int spinds_int
        row_start = row_end + 1;
    end
    A((J.f+1):end,:) = [Asl,Ass];
    A = sparse(A);
    c = sparse(c);
end