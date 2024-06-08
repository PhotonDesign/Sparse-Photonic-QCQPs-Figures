function[dDdn] = dDMatrix(lp, sp, sz)
% lp = layerProperties() struct
% sz = [Nw Nth]
I = ones(sz);
Z = zeros(sz);

dDdn = zeros([sz 4]);
if (strcmp(sp,'s'))
    dDdn = StampRow3(dDdn, Z, Z, lp.dkxdn, -lp.dkxdn);
else
    d_kn = lp.dkxdn./lp.n - lp.kx./lp.n.^2;
    dDdn = StampRow3(dDdn, d_kn, d_kn, I, -I);
end
end

