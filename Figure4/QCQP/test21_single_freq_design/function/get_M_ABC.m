function M = get_M_ABC(x, w)
    % actual M(x) divided by a factor of k^2
    Nx = length(x);
    dx = x(2) - x(1);
    M = sparse(1:Nx,1:Nx,-2*ones(Nx,1),Nx,Nx) + ...
        sparse(2:Nx,1:Nx-1,ones(Nx-1,1),Nx,Nx) + ...
        sparse(1:Nx-1,2:Nx,ones(Nx-1,1),Nx,Nx);
    M = 1/dx^2/w^2 * M;
    M = M + sparse(1:Nx,1:Nx,ones(Nx,1),Nx,Nx);
    M(1,:) = 0;
    M(1,[1,2]) = 1/dx^2/w^2 * [1j*w*dx-1, 1];
    M(Nx,:) = 0;
    M(Nx,[Nx-1,Nx]) = 1/dx^2/w^2 * [1, 1j*w*dx-1];
end