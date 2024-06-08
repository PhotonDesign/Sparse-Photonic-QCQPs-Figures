function g = GreensFunction_Kernel(k, ax, ay, hx, hy, quad_order)
    nx = ax / hx + 1;
    ny = ay / hy + 1;
    
    x2 = -ax:hx:ax; 
    y2 = -ay:hy:ay;
    [X2, Y2] = meshgrid(x2,y2);
    r = sqrt(X2.^2 + Y2.^2);
    g = 1j*k^2 / 4 * besselh(0,k*r);
    g(ny,nx) = 0; % set singularity (at the origin) to zero

    % coefficients for high-order quadrature correction
    if hx ~= hy 
        warning('Current high-order quadrature correction requires: ax / (nx-1) = ay / (ny-1)')
    end

    % high-order quadrature correction
    [D0, D1] = quad_D(k, hx);
    Tao = zeros(size(g));
    switch quad_order
        case 0 % no correction, do nothing
        case 4 % forth order correction
            Tao(ny,nx) = D0;
        case 6 % six order correction
            Tao(ny,nx) = D0 - 2*D1/hx^2;
            Tao(ny-1,nx) = 1/2 * D1 / hx^2;
            Tao(ny+1,nx) = 1/2 * D1 / hx^2;
            Tao(ny,nx-1) = 1/2 * D1 / hx^2;
            Tao(ny,nx+1) = 1/2 * D1 / hx^2;
    end
    g = g + 1j*k^2/4/hx^2 * Tao;
end

function [D0, D1] = quad_D(k, hx)
    [D0_t, D1_t] = planck_win_get_D(k, hx);
    D0 = D0_t * hx^2; 
    D1 = D1_t * hx^4;
end
