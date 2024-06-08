function [PML_bool,PML_bool_vect,PMLx_num,PMLy_num] = get_PML_bool_v2(Kx,Ky,K_PML)

PMLx_num = K_PML;%round(thickness*Kx/Lx);
PMLy_num = K_PML;%round(thickness*Ky/Ly);

PML_bool = zeros(Ky,Kx);
PML_bool(1:PMLy_num,:) = 1;
PML_bool(end-PMLy_num+1:end,:) = 1;
PML_bool(:,1:PMLx_num) = 1;
PML_bool(:,end-PMLx_num+1:end) = 1;
PML_bool_vect = reshape(PML_bool.',[],1);