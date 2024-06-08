
% This code computes the reflection and transmission properties of a 
% film with arbitrarily many layers.
% Premise: two half-spaces with N layers in between (thus: N+2 materials)
%   e.g. Air-SiO2-Si-SiO2-Air would be three layers, surrounded by two air
%   half-spaces
% NOTE: layer_thickness * frequency needs to be dimensionless
% Constraints: 
%    -- non-evanescent waves
%    -- initial/final layers assumed to be non-absorbing
% Implementation follows Yeh, Optical Waves in Layered Media, sec. 5.1

function[R,T,A,r,t,M] = multilayer_film(layerMaterials, layerThicknesses, w, theta, pol)

N = length(layerMaterials);
if (N ~= length(layerThicknesses)+2)
    error('Need to specify N layer thicknesses and N+2 materials');
end

Nw = length(w);
Nth = length(theta);
sz = [Nw Nth];
n0 = sqrt_k(layerMaterials{1}(w));
beta = repmat(n0(:),1,Nth) .* repmat(sin(theta(:).'),Nw,1); % normalized by c/w

R = cell(2,1); T = cell(2,1); A = cell(2,1);
r = cell(2,1); t = cell(2,1);
for p = 1:length(pol)
    M = MMatrix(layerMaterials, layerThicknesses, w, beta, pol(p), sz);
    
    M11 = M(:,:,1);
    M21 = M(:,:,3);
    r{p} = M21 ./ M11;
    t{p} = 1 ./ M11;

    lp0 = layerProps(sqrt_k(layerMaterials{1}(w)),beta,Nw,Inf); %kx0
    lps = layerProps(sqrt_k(layerMaterials{end}(w)),beta,Nw,Inf); %kxs
    R{p} = abs(r{p}).^2;
    T{p} = real(lps.kx)./lp0.kx .* abs(t{p}).^2;
    A{p} = 1 - R{p} - T{p};
end

end
