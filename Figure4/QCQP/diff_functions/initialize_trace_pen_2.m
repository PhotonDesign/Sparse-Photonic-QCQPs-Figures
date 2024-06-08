function [A,b,c,K,J] = initialize_trace_pen_2(F0,F,e)

    dim = length(F0);
    Nx = (dim-1)/2;
    
    f0 = reshape(F0,dim^2,1);
    f = cell(1,length(F));
    for i = 1:length(F)
        f{i} = reshape(F{i},length(F{i})^2,1);
    end
    
    b = sparse(length(f),1,1,length(f),1);
    
    tr_pen = e*reshape(eye(dim),dim^2,1);
    %Makes objective "min trace"
    
    c = sparse(f0 + tr_pen);
    
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
    
    A = sparse(row_idx_A, col_idx_A, val_A, length(f), dim^2);
    
    K.s = dim;

end