
% w, theta are scalars, not arrays
% AB is a 2*(N+2) x 1 array
% x is an Nx x 1 array
% f is an Nx x 1 array
function[f] = multilayer_fields_from_coeffs(x, AB, layerMaterials, layerThicknesses, w, theta)

n0 = sqrt_k(layerMaterials{1}(w));
beta = n0 * sin(theta); % normalized by c/w

Nw = length(w);
N = length(layerThicknesses);
xi = [0; cumsum(layerThicknesses(:))]; % N + 1 interfaces
lt = [Inf; layerThicknesses; Inf];

x = x(:);
f = zeros(length(x),1);
i = 1;
lp0 = layerProps(sqrt_k(layerMaterials{i}(w)), beta, Nw, lt(i)*w);
f = f + (AB(2*i-1)*exp(1i*lp0.kx*w*(x-xi(i))) + AB(2*i)*exp(-1i*lp0.kx*w*(x-xi(i)))).*(x<xi(i));
for i = 2:(N+1)
    lp = layerProps(sqrt_k(layerMaterials{i}(w)), beta, Nw, lt(i)*w);
    f = f + (AB(2*i-1)*exp(1i*lp.kx*w*(x-xi(i))) + AB(2*i)*exp(-1i*lp.kx*w*(x-xi(i)))).*(x<xi(i)).*(x>xi(i-1));
end
i = N+2;
lpS = layerProps(sqrt_k(layerMaterials{i}(w)), beta, Nw, lt(i)*w);
f = f + (AB(2*i-1)*exp(1i*lpS.kx*w*(x-xi(i-1))) + AB(2*i)*exp(-1i*lpS.kx*w*(x-xi(i-1)))).*(x>xi(i-1));
end
