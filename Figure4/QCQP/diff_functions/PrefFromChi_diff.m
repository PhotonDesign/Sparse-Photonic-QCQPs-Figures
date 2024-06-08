function [Pref,phase] = PrefFromChi_diff(chi,L,Nx)

    w = 1;
    L = L*2*pi;
    obj = 'R'; % option: R, Pext
    
    %%% init
    k = w;

    % material
    xi = - 1 ./ chi;

    % coordinate
    dx = L / (Nx-4);
    x = linspace(-2*dx, L+dx, Nx).';
    dx = x(2) - x(1);

    % Maxwell operator (with absorbing boundary conditions)
    M = get_M_ABC(x, w);

    % source
    rho = sparse(2, 1, 1/dx, Nx, 1);

    % calculate incident field 
    einc = M \ (- rho);
    
    e1 = (M + diag(chi)) \ (- rho);
    escat1 = e1 - einc;
    r = escat1(3) / einc(3); 
    Pref = abs(r)^2;
    
    phase = angle(einc(3)) - angle(escat1(3));
    if phase > pi
        phase = phase - 2*pi;
    elseif phase < -pi
        phase = phase + 2*pi;
    end
    
end