function G = GreensMatrix(k, xy)
    % Construct 2D, TE Green's function matrix without quadrature correction
    
    % Hankel function matrix
    Nxy = length(xy);
    H = zeros(Nxy,Nxy);
    for i = 1:Nxy
        dr = xy(i,:) - xy;
        dr = sqrt(dr(:,1).^2 + dr(:,2).^2);
        H(i,:) = besselh(0,dr);
    end

    Ngrid = length(xy);
    H(1:1+Ngrid:end) = 0; % set singularity to zero
    G = 1j*k^2 / 4 * H; % Green's functon matrix
end