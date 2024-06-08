
% Should just post-multiply by 0.5
function[dDIdn] = dDInvMatrix(lp, sp, sz)
Z = zeros(sz);
dDIdn = zeros([sz 4]);
if (strcmp(sp,'s'))
    dDIdn = StampRow3(dDIdn, Z, -0.5*lp.dkxdn./lp.kx.^2, Z, 0.5*lp.dkxdn./lp.kx.^2);
else
    d_nk = 1./lp.kx - lp.n./lp.kx.^2.*lp.dkxdn;
    dDIdn = StampRow3(dDIdn, 0.5*d_nk, -0.5./lp.n.^2, 0.5*d_nk, 0.5./lp.n.^2);
end
end