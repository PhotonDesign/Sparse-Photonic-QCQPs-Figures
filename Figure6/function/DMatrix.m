function[D] = DMatrix(lp, sp, sz)
% lp = layerProperties() struct
% sz = [Nw Nth]
I = ones(sz);
Z = zeros(sz);

D = zeros([sz 4]);
if (strcmp(sp,'s'))
    D = StampRow3(D, I, I, lp.kx, -lp.kx);
else
    D = StampRow3(D, lp.kx./lp.n, lp.kx./lp.n, lp.n, -lp.n);
end
end

