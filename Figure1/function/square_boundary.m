function xy_square = square_boundary(corner, Nx, Ny)    
    xy_square = [linspace_no_end(corner(1,1), corner(2,1), Nx-1)', corner(1, 2) * ones(Nx-1, 1);
                corner(2, 1) * ones(Ny-1, 1), linspace_no_end(corner(2,2), corner(3,2), Ny-1)';
                linspace_no_end(corner(3,1), corner(4,1), Nx-1)', corner(3, 2) * ones(Nx-1, 1);
                corner(4, 1) * ones(Ny-1, 1), linspace_no_end(corner(4,2), corner(1,2), Ny-1)';];
end


function y = linspace_no_end(a, b, n)
    y = linspace(a, b, n+1);
    y = y(1:n);
end