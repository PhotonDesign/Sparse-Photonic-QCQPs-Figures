
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

function[r,t,AB] = multilayer_film_coeffs(layerMaterials, layerThicknesses, w, theta, pol)

N = length(layerMaterials) - 2;
if (N ~= length(layerThicknesses))
    error('Need to specify N layer thicknesses and N+2 materials');
end

Nw = length(w);
Nth = length(theta);
sz = [Nw Nth];
n0 = sqrt_k(layerMaterials{1}(w));
beta = repmat(n0(:),1,Nth) .* repmat(sin(theta(:).'),Nw,1); % normalized by c/w

r = cell(2,1); t = cell(2,1);
AB = cell(2,1);
for p = 1:length(pol)
    M = MMatrix(layerMaterials, layerThicknesses, w, beta, pol(p), sz);
    
    M11 = M(:,:,1);
    M21 = M(:,:,3);
    r{p} = M21 ./ M11;
    t{p} = 1 ./ M11;
    
    ABp = zeros(Nw,Nth,2*(N+2));
    ABp(:,:,end-1) = t{p};
    ABp(:,:,end)   = 0;
    
    lt = [Inf; layerThicknesses; Inf];
    % interface N+1
    i = N+1;
    lp = layerProps(sqrt_k(layerMaterials{i}(w)), beta, Nw, lt(i)*w);
    lpp1 = layerProps(sqrt_k(layerMaterials{i+1}(w)), beta, Nw, lt(i+1)*w);
    DI = DInvMatrix(lp, pol, sz);
    D = DMatrix(lpp1, pol, sz);
    ind = (1:2) + 2*N;
    ABp(:,:,ind) = MVTimes(MTimes(DI,D), ABp(:,:,ind+2));
    % interface 1 through N
    for i = N:-1:1
        lp   = layerProps(sqrt_k(layerMaterials{i}(w)),   beta, Nw, lt(i)*w);
        lpp1 = layerProps(sqrt_k(layerMaterials{i+1}(w)), beta, Nw, lt(i+1)*w);
        DinvDP = DinvDPProd(lp, lpp1, pol, sz);
        ind = (1:2) + 2*(i-1);
        ABp(:,:,ind) = MVTimes(DinvDP, ABp(:,:,ind+2));
    end
    AB{p} = ABp;
end

end
