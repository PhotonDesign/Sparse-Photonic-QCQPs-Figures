function [Psca,RoverT] = getPscaFromChi(chi,L,M)
    c = 1;
    wl = 1;
    f = c/wl;
    w = 2*pi*f;
    % Units normalized to 1
    zMat = linspace(0,L,M)';
    % Pick desired length and discretization
    
    EMat = get_E(zMat,2*pi);
    chiMat = diag(1./chi);
    GMat = zeros(M);
    for a = 1:M
        for b = 1:M
            GMat(a,b) = get_G(zMat(a),zMat(b),2*pi,L/M);
        end
    end
    
    GMatExt2 = zeros(1,M); %Higher resolution grid to account for phase
    for b = 1:M
        GMatExt2(b) = get_G(-0.5,zMat(b),2*pi,L/M); %arbitrary point behind scatterer
    end
    GMatExt1 = zeros(1,M);
    for b = 1:M
        GMatExt1(b) = get_G(L+0.5,zMat(b),2*pi,L/M); %arbitrary point ahead of scatterer
    end
    
    phiMat = linsolve(GMat-chiMat,-EMat);
    Pref = abs(GMatExt2*phiMat)^2;
    Pfwd = abs(GMatExt1*phiMat)^2;
    Psca = Pref + Pfwd;
    RoverT = Pref/Pfwd;
end