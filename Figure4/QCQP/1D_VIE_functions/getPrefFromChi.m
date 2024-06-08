function [Pref,Eref,phiMat,phase] = getPrefFromChi(chi,L,M,wOverw0)
    c = 1;
    w = 2*pi*wOverw0;
    k = w/c;
    % Units normalized to 1
    zMat = linspace(0,L,M)';
    %L is given in units of c/w0 (w0 vacuum wavelength)
    % Pick desired length and discretization
    
    EMat = get_E(zMat,k);
    chiMat = diag(1./chi);
    GMat = zeros(M);
    for a = 1:M
        for b = 1:M
            GMat(a,b) = get_G(zMat(a),zMat(b),k,L/M);
        end
    end
    
    N = 1000;
    zExt = linspace(-1,0,N);
    GMatExt2 = zeros(N,M); %Higher resolution grid to account for phase
    for a = 1:N
        for b = 1:M
            GMatExt2(a,b) = get_G(zExt(a),zMat(b),k,L/M); %arbitrary point behind scatterer
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
end