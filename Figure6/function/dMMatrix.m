function[dM] = dMMatrix(m, t, w, beta, pol, sz)
% m = layerMaterials, t = layerThicknesses
% sz = [Nw Nth]
% This function computes dMdn, *not* dMdEps
Nw = sz(1);
N = length(m);

% exit half-space
lp = layerProps(sqrt(m{end}(w)), beta, Nw, Inf);
D = DMatrix(lp, pol, sz);
dM = repmat(D, [1 1 1 N-2]);

% interior layers
for i = N-1:-1:2
    lp = layerProps(sqrt(m{i}(w)), beta, Nw, t(i-1)*w);
    [DPD,dDPD] = DPDinvProd(lp, pol, sz, 1);
    dM(:,:,:,i-1) = MTimes(dDPD,dM(:,:,:,i-1));
    % oddly, this seems slightly faster than vectorized, below
    for j = [2:i-1 i+1:N-1] % j ~= i
        dM(:,:,:,j-1) = MTimes(DPD,dM(:,:,:,j-1));
    end
%     ind = [1:i-2 i:N-2];
%     dM(:,:,:,ind) = MTimes(repmat(DPD, [1 1 1 N-3]), dM(:,:,:,ind));
end

% entrance half-space
lp = layerProps(sqrt(m{1}(w)),beta,Nw,Inf);
DInv = repmat(DInvMatrix(lp, pol, sz), [1 1 1 N-2]);
dM = MTimes(DInv, dM);
end