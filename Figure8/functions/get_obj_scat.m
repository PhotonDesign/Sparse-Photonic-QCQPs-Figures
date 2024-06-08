function [S0,b0,c0] = get_obj_scat(obj,Wm,MM_bg,omega,k,eps0,einc_vect)

%     KxKy = size(einc_vect,1);
      
if ~strcmp(obj,'scattering')
    disp('Wrong obj.')
    return;
    
else   
    S0 = Im((omega/2)*(eps0/k^2)*Wm*MM_bg);
    b0 = 0*einc_vect;
    c0 = 0;
    
    S0 = sparse(S0);
    b0 = sparse(b0);
    c0 = sparse(c0);
end