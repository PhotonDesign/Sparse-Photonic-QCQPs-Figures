function [S0,b0,c0] = get_obj_e_cent(obj,center_idx,einc_vect)
    
    KxKy = size(einc_vect,1);
    
    center_mat = zeros(KxKy,KxKy);
    center_mat(center_idx,center_idx) = 1;
   
if ~strcmp(obj,'full field at center')
    disp('Wrong obj.')
    return;
    
else   
    S0 = center_mat;
    b0 = 2*center_mat'*einc_vect;
    c0 = abs(einc_vect(center_idx))^2;
    
    S0 = sparse(S0);
    b0 = sparse(b0);
    c0 = sparse(c0);
end