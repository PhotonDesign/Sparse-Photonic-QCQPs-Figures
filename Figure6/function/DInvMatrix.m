
% Should just post-multiply by 0.5
function[Dinv] = DInvMatrix(lp, sp, sz)
I = ones(sz);
Dinv = zeros([sz 4]);
if (strcmp(sp,'s'))
    Dinv = StampRow3(Dinv, 0.5*I, 0.5./lp.kx, 0.5*I, -0.5./lp.kx);
else
    Dinv = StampRow3(Dinv, 0.5*lp.n./lp.kx, 0.5./lp.n, 0.5*lp.n./lp.kx, -0.5./lp.n);
end
end