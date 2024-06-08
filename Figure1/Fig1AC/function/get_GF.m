function G = get_GF(k, xy_out, xy_in)
    Gr = @(r) 1j*k^2 / 4 * besselh(0,k*r);
    Nout = size(xy_out,1);
    Nin = size(xy_in,1);
    G = zeros(Nout, Nin);
    for i = 1:Nout
        dr = vecnorm(xy_in - xy_out(i, :), 2, 2);
        G(i, :) = Gr(dr);
    end
end