function[A_VIE,b_VIE,c_VIE] = get_obj_focus_VIE2(obj,G_tar_VIE,einc_tar)

if ~strcmp(obj,'focusing')
    disp('Wrong obj.')
    return;
    
else
    
KxKy_VIE = size(G_tar_VIE,2);
             
%% get VIE target function

% A_VIE = zeros(KxKy_VIE,KxKy_VIE);
% b_VIE = exp(-1i*theta)*G_tar_VIE';
% c_VIE = real(einc_tar*exp(1i*theta));  


A_VIE = G_tar_VIE'*G_tar_VIE;

b_VIE = 2*G_tar_VIE'*einc_tar;

c_VIE = abs(einc_tar)^2;


end