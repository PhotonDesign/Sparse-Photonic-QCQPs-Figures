
% obj = 'R'; % option: R, Pext

% targetPhase = 0.3*pi;
% phase_dif = pi/2 - targetPhase; %Differential starts out of phase w.r.t VIE


function[etot,escat,einc,r,R] = fdfd_1d(x,dx,chi,k,Psrc_amp)

Nx = length(x);
M = get_M_ABC(x, k); % Maxwell operator with ABC

if (nargin == 4)
    Psrc_amp = 2/(k*dx); % plane-wave with near-unity amplitude (normalized out for reflection coefficient anyhow)
end
Psrc_amp = Psrc_amp * exp(-1i*k*dx);
rho = sparse(2, 1, Psrc_amp, Nx, 1);

einc = M \ (- rho);
etot = (M + diag(chi)) \ (- rho);
escat = etot - einc;

r = escat(3) / einc(3); 
R = abs(r)^2;
end

