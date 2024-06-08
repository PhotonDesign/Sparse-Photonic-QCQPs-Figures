function e = film_current(epsr,z)
    % calculate analytical solution of current in a 1D film
    % film material: epsr
    % film coordinate: z
    % assume: normal unit-amplitude plane wave incidence
    n = sqrt(epsr);
    z1 = z(1);
    z2 = z(end);
    A2 = 2*(n+1)*exp(1j*z1)*exp(1j*n*(z1-2*z2)) ...
        / (-(n-1)^2 + (n+1)^2*exp(1j*n*2*(z1-z2)));
    B2 = 2*(n-1)*exp(1j*z1)*exp(1j*n*z1) ...
        / (-(n-1)^2 + (n+1)^2*exp(1j*n*2*(z1-z2)));
    e = A2*exp(1j*n*z) + B2*exp(-1j*n*z);
end