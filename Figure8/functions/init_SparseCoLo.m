function[S, b, c] = init_SparseCoLo(bg_idx_vect,des_idx_vect,MM_bg,xi,k,einc_vect)

KxKy = length(bg_idx_vect) + length(des_idx_vect);

M0 = sparse(KxKy,KxKy);
M0(des_idx_vect,:) = MM_bg(des_idx_vect,:);

S = {};
b = {};
c = {};
Mquad = sparse(KxKy,KxKy);

for ii = bg_idx_vect.'    

    S(end+1) = {sparse(KxKy, KxKy)};
    b(end+1) = {MM_bg(ii,:)'};
    c(end+1) = {sparse(0)};
    
    Mquad = Mquad + MM_bg(ii,:)' * MM_bg(ii,:); % summed redundant quadratic constraints (background region)
end

for ii = des_idx_vect.'
    
    Di = sparse(KxKy, KxKy); % D-matrix for points in the degsign region
    Di(ii,ii) = 1;
    
    S(end+1) = {Di * M0 - (xi'./k^2) * (M0' * Di * M0)};       
    b(end+1) = {M0' * Di * einc_vect};
    c(end+1) = {sparse(0)};
end

% summed redundant quadratic constraints (background region)
S(end+1) = {Mquad};
b(end+1) = {sparse(KxKy,1)};
c(end+1) = {sparse(0)};

% complex to real constraints
[S, b, c] = c2r_con(S, b, c);