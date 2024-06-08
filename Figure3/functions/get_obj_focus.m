function [S0,b0,c0] = get_obj_focus(obj,theta,G_tar,MM_bg,omega,einc_tar)


if ~strcmp(obj,'focusing')
    disp('Wrong obj.')
    return;
    
else
    
    KxKy = size(MM_bg,1);
    
    S0 = zeros(KxKy,KxKy);
    b0 = exp(-1i*theta)*(-1/omega^2)*MM_bg'*G_tar';
    c0 = real(einc_tar*exp(1i*theta));
    
    S0 = sparse(S0);
    b0 = sparse(b0);
    c0 = sparse(c0);
end