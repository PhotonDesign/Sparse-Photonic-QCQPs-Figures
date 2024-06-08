% Created on Feb 11, 2019 by Wenjin Xue

function[dM] = dMMatrix_thickness(m, t, w, beta, pol, sz)
% m = layerMaterials, t = layerThicknesses
% sz = [Nw Nth]
% This function computes dMdt
Nw = sz(1);
N = length(m);

% exit half-space
lp = layerProps(sqrt(m{end}(w)), beta, Nw, Inf);
D = DMatrix(lp, pol, sz);
dM = repmat(D, [1 1 1 N-2]);

% interior layers
for i = N-1:-1:2
    lp = layerProps(sqrt(m{i}(w)), beta, Nw, t(i-1)*w);
    [DPD,dDPD] = DPDinvProd_thickness(lp, pol, sz,w, 1);
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

% dFlag = 1 if want deriv's
function[DPD, dDPDdt] = DPDinvProd_thickness(lp, pol, sz,w, dFlag)
D = DMatrix(lp, pol, sz);
DI = DInvMatrix(lp, pol, sz);
P = PMatrix(lp, sz);
DPD = Times3(D,P,DI);

if (dFlag)
    Z = zeros(sz);
    dPdt = zeros([sz 4]);
    dexpdt = -1i*w*lp.kx.*exp(-1i*lp.kx.*lp.t);
    dPdt = StampRow3(dPdt, dexpdt, Z, Z, conj(dexpdt));
    dDPDdt = Times3(D,dPdt,DI);
end
end