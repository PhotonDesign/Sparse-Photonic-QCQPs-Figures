function[S, b, c] = init_SparseCoLo_giveD(bg_idx_vect,des_idx_vect,MM,xi,k,einc_vect,D)

KxKy = length(bg_idx_vect) + length(des_idx_vect);

M0des = sparse(KxKy,KxKy);
M0des(des_idx_vect,:) = MM(des_idx_vect,:);

M0bg = sparse(KxKy,KxKy);
M0bg(bg_idx_vect,:) = MM(bg_idx_vect,:);

S = {};
b = {};
c = {};


% background quadratic
S(end+1) = {M0bg' * M0bg};
b(end+1) = {sparse(KxKy,1)};
c(end+1) = {sparse(0)};

for tt = 1:length(D)
    
Dtt = sparse(D{tt});

Ddes = sparse(KxKy,KxKy);
Ddes(des_idx_vect,:) = Dtt(des_idx_vect,:);

Dbg = sparse(KxKy,KxKy);
Dbg(bg_idx_vect,:) = Dtt(bg_idx_vect,:);

% background  
S(end+1) = {sparse(KxKy,KxKy)};
b(end+1) = {( ones(KxKy,1)' * Dbg * M0bg )'};
c(end+1) = {sparse(0)};

% design region
S(end+1) = {Ddes * M0des - (xi'./k^2) * (M0des' * Ddes * M0des)};
b(end+1) = {( einc_vect' * Ddes * M0des )'};
c(end+1) = {sparse(0)};


% tmp_idx = find(diag(Dtt));
% if find(bg_idx_vect==tmp_idx) % background  
%   
% S(end+1) = {sparse(KxKy,KxKy)};
% b(end+1) = {( ones(KxKy,1)' * Dtt * M0bg )'};
% c(end+1) = {sparse(0)};
% 
% else % design region
% S(end+1) = {Dtt * M0des - (xi'./k^2) * (M0des' * Dtt * M0des)};
% b(end+1) = {( einc_vect' * Dtt * M0des )'};
% c(end+1) = {sparse(0)};
% 
% end

end

[S, b, c] = c2r_con(S, b, c);