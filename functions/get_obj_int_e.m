function [S0,b0,c0] = get_obj_int_e(obj,Wm,einc_vect)

    KxKy = size(einc_vect,1);
   
if ~strcmp(obj,'integrate full field') && ~strcmp(obj, 'full field at center line')...
        && ~strcmp(obj,'integrate full field rect') 
    disp('Wrong obj.')
    return;
    
else   
    S0 = eye(KxKy,KxKy)*Wm;
    b0 = 2*Wm'*einc_vect;
    c0 = einc_vect'*Wm*einc_vect;
    
    S0 = sparse(S0);
    b0 = sparse(b0);
    c0 = sparse(c0);
end