function Wm = get_trapz_mat(Kx,Ky,Nx,Mx,Ny,My)

% 2D trapezoidal integration matrix for design region

W = zeros(Ky,Kx);

W(Ny+1:My,Nx+1:Mx) = 1;
W(Ny+1,Nx+1:Mx) = 1/2;
W(Ny+1:My,Nx+1) = 1/2;
W(My,Nx+1:Mx) = 1/2;
W(Ny+1:My,Mx) = 1/2;
W(Ny+1,Nx+1) = 1/4;
W(Ny+1,Mx) = 1/4;
W(My,Nx+1) = 1/4;
W(My,Mx) = 1/4;

W_vect = reshape(W.',[],1);

Wm = diag(W_vect);