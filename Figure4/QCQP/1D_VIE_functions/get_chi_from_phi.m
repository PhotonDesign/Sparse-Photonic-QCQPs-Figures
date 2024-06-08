function chiMat = get_chi_from_phi(L,M,phiMat)
    zMat = linspace(0,L,M)';
    Einc = get_E(zMat,2*pi);
    GMat = zeros(M,M);
    for a = 1:M
        for b = 1:M
            GMat(a,b) = get_G(zMat(a),zMat(b),2*pi,L/M);
        end
    end
    A = diag(Einc + GMat*phiMat);
    chiMat = linsolve(A,phiMat);
end