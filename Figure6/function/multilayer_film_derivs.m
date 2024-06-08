
% This code computes *derivatives* with respect to the reflection and 
% transmission properties of a film with arbitrarily many layers.
% Premise: two half-spaces with N layers in between (thus: N+2 materials)
%   e.g. Air-SiO2-Si-SiO2-Air would be three layers, surrounded by two air
%   half-spaces
% NOTE: layer_thickness * frequency needs to be dimensionless
% Constraints: 
%    -- non-evanescent waves
%    -- initial/final layers assumed to be non-absorbing
% Forward solver implementation follows Yeh, Optical Waves in Layered Media, sec. 5.1
%
% Outputs: size Nw x Nt x 1 x N, where N is # of refractive indices

function[R,dR,r,dr,T,dT,t,dt] = multilayer_film_derivs(layerMaterials, layerThicknesses, w, theta, pol)

N = length(layerMaterials);
if (N ~= length(layerThicknesses)+2)
    error('Need to specify N layer thicknesses and N+2 materials');
end

Nw = length(w);
Nth = length(theta);
sz = [Nw Nth];
n0 = sqrt_k(layerMaterials{1}(w));
beta = repmat(n0(:),1,Nth) .* repmat(sin(theta(:).'),Nw,1); % normalized by c/w

R = cell(2,1); T = cell(2,1);
dR = cell(2,1); dT = cell(2,1);
r = cell(2,1); t = cell(2,1);
dt = cell(2,1); dr = cell(2,1);
for p = 1:length(pol)
    M = MMatrix(layerMaterials, layerThicknesses, w, beta, pol(p), sz);
    M11 = repmat(M(:,:,1), [1 1 1 N-2]);
    M21 = repmat(M(:,:,3), [1 1 1 N-2]);

    dMdn = dMMatrix(layerMaterials, layerThicknesses, w, beta, pol(p), sz);
    dM11 = dMdn(:,:,1,:);
    dM21 = dMdn(:,:,3,:);
    
    r{p} = M21 ./ M11;
    R{p} = abs(r{p}).^2;
    dR{p} = 2./abs(M11).^2 .* real(dM21.*conj(M21) - R{p}.*dM11.*conj(M11));
    dr{p} = dM21 ./ M11 - r{p} .* dM11 ./ M11;
    R{p} = R{p}(:,:,1,1);
    r{p} = r{p}(:,:,1,1);
    
    if (nargout>4)
        lp0 = layerProps(sqrt(layerMaterials{1}(w)),beta,Nw,Inf); %kx0
        lps = layerProps(sqrt(layerMaterials{end}(w)),beta,Nw,Inf); %kxs
        kx0 = repmat(lp0.kx, [1 1 1 N-2]);
        kxs = repmat(lps.kx, [1 1 1 N-2]);
        t{p} = sqrt(kxs./kx0) ./ M11;
        T{p} = abs(t{p}).^2;
        %T{p} = kxs./kx0 .* abs(1./M11).^2;
        dT{p} = -2*real(dM11./M11) .* T{p};
        dt{p} = -t{p} .* dM11./M11;
        T{p} = T{p}(:,:,1,1);
        t{p} = t{p}(:,:,1,1);
    end

%     if (nargout>4)
%         A{p} = 1 - R{p} - T{p};
%         dA{p} = -dR{p} - dT{p};
%     end
end

end
