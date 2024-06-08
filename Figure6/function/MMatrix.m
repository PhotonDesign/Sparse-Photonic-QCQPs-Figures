function[M] = MMatrix(m, t, w, beta, pol, sz)
% m = layerMaterials, t = layerThicknesses
% sz = [Nw Nth]
Nw = sz(1);
Nm = length(m);

% exit half-space
lp = layerProps(sqrt_k(m{end}(w)), beta, Nw, Inf);
M = DMatrix(lp, pol, sz);

% interior layers
for i = Nm-1:-1:2
    lp = layerProps(sqrt_k(m{i}(w)), beta, Nw, t(i-1)*w);
    M = MTimes(DPDinvProd(lp, pol, sz, 0), M);
end

% entrance half-space
lp = layerProps(sqrt_k(m{1}(w)),beta,Nw,Inf);
M = MTimes(DInvMatrix(lp, pol, sz), M);
end