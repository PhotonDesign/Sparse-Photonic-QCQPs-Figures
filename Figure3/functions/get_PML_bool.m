function [PML_bool,PML_bool_vect,PMLx_num,PMLy_num] = get_PML_bool(thickness,Kx,Ky,Lx,Ly)

PMLx_num = round(thickness*Kx/Lx);
PMLy_num = round(thickness*Ky/Ly);

PML_bool = zeros(Ky,Kx);
PML_bool(1:PMLy_num,:) = 1;
PML_bool(end-PMLy_num+1:end,:) = 1;
PML_bool(:,1:PMLx_num) = 1;
PML_bool(:,end-PMLx_num+1:end) = 1;
PML_bool_vect = reshape(PML_bool.',[],1);