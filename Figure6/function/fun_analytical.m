
function [xx, field] = fun_analytical(Nx)



rng(13)
Nm = 10; % number of material patches
L = 8*pi+0.4;
xm = sort(rand(2*Nm,1)) * L;
xm(1) = 0.2;
xm(end) = xm(1)+8*pi;
xm(end-5) = 19.2;



%%
% geometry:
%   vacuum, const n, vacuum
%   Yeh, p. 91

d = diff(xm);
eps = 12;
const_epsr1 = @(w) eps*ones(size(w));
const_epsr2 = @(w) 1^2*ones(size(w));
layer_incide = {@Vacuum};
for i = 1:Nm-1
    layer_incide{end+1} = const_epsr1;
    layer_incide{end+1} = const_epsr2;
end
layer_incide{end+1} = const_epsr1;
layer_incide{end+1} = @Vacuum;
layerMaterials = layer_incide;

pol = {'s'}; % could also be {'s'} or {'p','s'}

w = 1;
theta = 0;
% d = [0.389; 0.2; 0.5];

[r,t,AB] = multilayer_film_coeffs(layerMaterials, d, w, theta, pol);

tot_thick = xm(end)-xm(1);
x = linspace(-tot_thick/(Nx-1),1.000*tot_thick + tot_thick/(Nx-1),Nx+2).';
f = multilayer_fields_from_coeffs(x, AB{1}, layerMaterials, d, w, theta);

%%

xx = x(2:end-2);
field = f(2:end-2);


end