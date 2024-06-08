function[Dbg,Ddes] = get_next_D_DIFF(es_opt,bg_idx_vect,des_idx_vect,M,xi,k,einc_vect)


es_opt = es_opt(:); % make sure it's a column
einc_vect = einc_vect(:); % make sure it's a column


KxKy = length(bg_idx_vect) + length(des_idx_vect);

M0des = sparse(KxKy,KxKy);
M0des(des_idx_vect,:) = M(des_idx_vect,:);

M0bg = sparse(KxKy,KxKy);
M0bg(bg_idx_vect,:) = M(bg_idx_vect,:);



% design
d_des = (M0des*(es_opt*es_opt') - (xi'./k^2)*M0des*(es_opt*es_opt')*M0des' + M0des*es_opt*einc_vect')';

Ddes = diag(diag(d_des));


% background
d_bg = ( M0bg * es_opt * ones(KxKy,1)' )';

Dbg = diag(d_bg);


% Ddes = Ddes./sqrt(sum(abs(diag(Ddes)).^2));
% 
% Dbg = Dbg./sqrt(sum(abs(diag(Dbg)).^2));



end