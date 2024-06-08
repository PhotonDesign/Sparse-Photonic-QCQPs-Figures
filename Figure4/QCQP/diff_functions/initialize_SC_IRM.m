function [A,b,c,K,J] = initialize_SC_IRM(F0,F,wk,X,bags,nrow_As)
    nb = length(bags);
    
    [V0,D0] = eigs(X);
    x0 = V0(:,1)*sqrt(D0(1,1));
    X0 = [X,x0;x0',1];
    Vk = submat_eigs_IRM_SC(X0,bags);
    
    F0 = sparse([F0,zeros(length(F0),1);zeros(1,length(F0)),0]);
    f0 = reshape(F0,length(F0)^2,1);
    f = cell(1,length(F)+1);
    for i = 1:length(F)
        F{i} = sparse([F{i},zeros(length(F{i}),1);zeros(1,length(F{i})),0]);
        f{i} = reshape(F{i},length(F{i})^2,1);
    end
    f{end} = sparse(length(f{1}),1,1,length(f{1}),1);
    
    oneInds = [length(f)-1,length(f)];
    b = [sparse(oneInds,1,1,length(f),1);sparse(nrow_As,1)];
    
    c = sparse([ones(nb,1)*wk;f0]);
    
    J.f = length(f);
    %J.s = repelem(tw+1,nb);
    
    %K.l = nb;
    K.s = [ones(1,nb) length(F0)]; %F0 has same dims as X variable
                      %x is same length as f0 (length(F0)^2)
    
    %Make Afs
    row_idx_A = [];
    col_idx_A = [];
    val_A = [];
    for i = 1:length(f)
        [row_idx,col_idx,val] = find(f{i}.');
        row_idx = row_idx + (i-1);
        row_idx_A = [row_idx_A, row_idx];
        col_idx_A = [col_idx_A, col_idx];
        val_A = [val_A, val];
    end
    
    Afs = sparse(row_idx_A, col_idx_A, val_A, length(f), length(f0));
    Afl = sparse(length(f),nb);
    
    %Make As. row (Asl and Ass)
%     Asl = zeros(psd_size,nb);
%     Ass = zeros(psd_size,length(f0));
    Asl = sparse(nrow_As,nb);
    Ass = sparse(nrow_As,length(f0));
    J.s = [];
    row_start = 1;
    for k = 1:nb
        cur_bag = [bags{k}, length(F0)];
        bag_size = length(cur_bag);
        noVs = bag_size - 1;
        row_end = row_start + noVs^2 - 1;
        Asl(row_start:row_end,k) = reshape(speye(noVs),noVs^2,1);
        Ass_int = zeros(noVs^2,bag_size^2);
        spinds_int = zeros(bag_size^2,length(f0));
        for i = 1:noVs
            for j = 1:noVs
                row = noVs*(i-1)+j;
                Ass_int(row,:) = -reshape(Vk{k}(:,i)*Vk{k}(:,j)',bag_size^2,1);
            end
        end
        for i = 1:bag_size
            for j = 1:bag_size
                row = bag_size*(i-1)+j;
                ind = sub2ind(size(F0),cur_bag(j),cur_bag(i));
                spinds_int(row,:) = sparse(1,ind,1,1,length(f0));
            end
        end
        Ass(row_start:row_end,:) = Ass_int*spinds_int;
        %spinds_int gives PSD constraint
        clear Ass_int spinds_int
        J.s = [J.s,noVs];
        row_start = row_end + 1;
    end
    
    A = sparse([Afl,Afs;Asl,Ass]);
    
 
    
        
    
    
    