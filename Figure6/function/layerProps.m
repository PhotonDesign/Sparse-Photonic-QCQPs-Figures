function[lp] = layerProps(n,beta,Nw,t)
lp = struct; % lp = layer prop's

% ref index
nx_size_inc = Nw / size(n,1);
Nth = size(beta,2);
n = repmat(n,nx_size_inc,Nth);
lp.n = n;

% kx (normalized to w/c)
kx = sqrt(n.^2 - beta.^2);
lp.kx = kx;
lp.dkxdn = n./kx;

% thickness
if (~isinf(t))
    t = repmat(t,1,Nth);
end
lp.t = t;
end