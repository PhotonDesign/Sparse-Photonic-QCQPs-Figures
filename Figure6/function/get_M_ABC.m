function M = get_M_ABC(x, w)
    Nx = length(x);
    dx = x(2) - x(1);
    M = -2 * diag(ones(Nx,1)) + diag(ones(Nx-1,1), 1) + diag(ones(Nx-1,1), -1);
    M = 1/dx^2 * M;
    M = M + w^2 * diag(ones(Nx,1));
    M(1,:) = 0;
    M(1,[1,2]) = 1/dx^2 * [1j*w*dx-1, 1];
    M(Nx,:) = 0;
    M(Nx,[Nx-1,Nx]) = 1/dx^2 * [1, 1j*w*dx-1];
end