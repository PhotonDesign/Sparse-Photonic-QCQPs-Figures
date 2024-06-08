function[dPdn] = dPMatrix(lp, sz)
Z = zeros(sz);
dPdn = zeros([sz 4]);
dexpNeg = -1i*lp.dkxdn.*lp.t.*exp(-1i*lp.kx.*lp.t);
dPdn = StampRow3(dPdn, dexpNeg, Z, Z, conj(dexpNeg));
end