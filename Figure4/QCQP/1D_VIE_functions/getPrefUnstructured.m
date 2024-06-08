function [Pref,Eref,phiMat,phase] = getPrefUnstructured(n,L,M)
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
    
    N = 1000;
    zExt = linspace(-1,0,N);
    GMatExt2 = zeros(N,M); %Higher resolution grid to account for phase
    for a = 1:N
        for b = 1:M
            GMatExt2(a,b) = get_G(zExt(a),zMat(b),2*pi,L/M); %arbitrary point behind scatterer
        end
    end
    
    phiMat = linsolve(GMat-chiMat,-EMat);
    Eref = real(GMatExt2*phiMat);
    Pref = abs(GMatExt2(1,:)*phiMat)^2;
        [~,maxInd] = max(Eref);
    phase = -zExt(maxInd)*2*pi;
    if phase > pi
        phase = phase-2*pi;
    end