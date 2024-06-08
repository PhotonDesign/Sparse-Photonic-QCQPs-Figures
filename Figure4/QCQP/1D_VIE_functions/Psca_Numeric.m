function PscatN = Psca_Numeric(n,L,M)
    c = 1;
    wl = 1;
    f = c/wl;
    w = 2*pi*f;
    % Units normalized to 1
    chi = n^2 - 1;
    zMat = linspace(0,L,M)';
    % Pick desired length and discretization
    
    EMat = get_E(zMat,2*pi);
    chiMat = (1/chi)*eye(M);
    GMat = zeros(M,M);
    for a = 1:M
        for b = 1:M
            GMat(a,b) = get_G(zMat(a),zMat(b),2*pi,L/M);
        end
    end
    
    coefMat = -(GMat-chiMat);
    pMat = linsolve(coefMat,EMat);
    PscatN = (w/2)*pMat'*(((GMat-GMat')/(2i))*pMat)*(L/M)*2;
end