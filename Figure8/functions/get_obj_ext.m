function [S0,b0,c0] = get_obj_ext(obj,Wm,MM_bg,omega,k,eps0,einc_vect)

    KxKy = size(einc_vect,1);
      
if ~strcmp(obj,'extinction')
    disp('Wrong obj.')
    return;
    
else   
    S0 = sparse(KxKy,KxKy);
    b0 = (omega/2).*(-1i)*(eps0/k^2)*MM_bg'*Wm'*einc_vect;
    c0 = 0;
    
    S0 = sparse(S0);
    b0 = sparse(b0);
    c0 = sparse(c0);
end