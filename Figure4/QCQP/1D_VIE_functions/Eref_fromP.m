function [Eref] = Eref_fromP(p,L)
    M = length(p);
    zExt = linspace(-2,0,M);
    zMat = linspace(0,L,M);
    GExt = zeros(M);
    for i = 1:M
        for j = 1:M
            GExt(i,j) = get_G(zExt(i),zMat(j),2*pi,L/M);
        end
    end
    Eref = GExt*p;
end