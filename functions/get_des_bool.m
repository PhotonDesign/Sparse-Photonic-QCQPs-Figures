function [des_bool,des_bool_vect,bg_idx_vect,des_idx_vect] = get_des_bool(Kx,Ky,Nx,Mx,Ny,My)

des_bool = zeros(Ky,Kx);
des_bool(Ny+1:My,Nx+1:Mx) = 1;
des_bool_vect = reshape(des_bool.', [], 1);

bg_idx_vect = find(~des_bool_vect);
des_idx_vect = find(des_bool_vect);
