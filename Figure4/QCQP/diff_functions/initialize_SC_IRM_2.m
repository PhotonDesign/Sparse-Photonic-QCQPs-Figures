function [A,b,c,K,J] = initialize_SC_IRM_2(F0,F,wk,X)
    
    [V0,D0] = eig(X);
    D0 = diag(D0);
    dim = length(D0);
    Vk = V0(1:dim,1:(dim-1));
    
    f0 = reshape(F0,length(F0)^2,1);
    f = cell(1,length(F));
    for i = 1:length(F)
        f{i} = reshape(F{i},length(F{i})^2,1);
    end
    
    b = sparse(length(f),1,1,length(f),1);
    %Increase size later
    
    c = sparse([wk;f0]);
    
    J.f = length(f);
    
    %K.l = nb;
    K.s = [1, length(F0)]; %F0 has same dims as X variable
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
    Afl = sparse(length(f),1);
    
    %Make As. row (Asl and Ass)
    [ncol,nrow] = size(Vk); %Relationship flips for A
    Ass_1 = reshape(speye(nrow),nrow^2,1);
    Ass = zeros(nrow^2,ncol^2);

    for i = 1:nrow
        for j = 1:nrow
            row = nrow*(i-1)+j;
            Ass(row,:) = -reshape(Vk(:,i)*Vk(:,j)',1,ncol^2);
        end
    end
    
    J.s = nrow;
    b = [b;sparse(nrow^2,1)];
    
    A = sparse([Afl,Afs;Ass_1,Ass]);
    b = sparse(b);
    
 
    
        
    
    
    