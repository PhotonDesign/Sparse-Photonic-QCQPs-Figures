
function[x,dx,Nx,designInd,Ndesign] = set_fdfd_grid(L,Nx,xDesignStart)
dx = L / (Nx-4);
x = linspace(-2*dx, L+dx, Nx).';
dx = x(2) - x(1);

[~,startInd] = min(abs(x-xDesignStart));
designInd = (startInd:Nx).';
Ndesign = length(designInd);
end