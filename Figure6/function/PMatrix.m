function[P] = PMatrix(lp, sz)
Z = zeros(sz);
P = zeros([sz 4]);
P = StampRow3(P, exp(-1i*lp.kx.*lp.t), Z, Z, exp(1i*lp.kx.*lp.t));
end