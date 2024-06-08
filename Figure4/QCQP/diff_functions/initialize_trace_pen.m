function [A,b,c,K,J] = initialize_trace_pen(F0,F,t)
    
    dim = length(F0);
    Nx = (dim-1)/2;

    f = cell(1,length(F));
    for i = 1:length(F)
        f{i} = reshape(F{i},length(F{i})^2,1);
    end
    
    bf = sparse(length(f),1,1,length(f),1);
    
    c = sparse(reshape(F0,dim^2,1));
    %Makes objective "min trace"
    
    J.f = length(F);
    
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
    
    Afs = sparse(row_idx_A, col_idx_A, val_A, length(f), dim^2);
    
    Ass = -reshape(speye(dim),1,dim^2);
    
    bs = -t;
    
    J.s = 1; %Scalar PSD cone
    
    K.s = dim;
    
    A = sparse([Afs;Ass]);
    b = sparse([bf;bs]);

end