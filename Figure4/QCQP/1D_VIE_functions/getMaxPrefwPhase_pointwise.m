function [maxP,phiMat,Eref,phase,solnMat] = getMaxPrefwPhase_pointwise(n,L,M,inputPhase)
    if (inputPhase < -pi || inputPhase > pi)
        error('Input phase must be between -pi and pi.') 
    end

    c = 1;
    wl = 1;
    f = c/wl;
    w = 2*pi*f;
    % Units normalized to 1
    chi = n^2 - 1;
    zMat = linspace(0,L,M)';
    % Pick desired length and discretization
    
    EMat = get_E(zMat,2*pi);
    chiInv = (1/chi)*eye(M);
    GMat = zeros(M,M);
    for a = 1:M
        for b = 1:M
            GMat(a,b) = get_G(zMat(a),zMat(b),2*pi,L/M);
        end
    end
   
    DMat = zeros(M,M,2*M+2);
    AMat = zeros(M,M,2*M+2);
    BMat = zeros(M,1,2*M+2);
    BMatT = zeros(1,M,2*M+2);
    COpt = zeros(M+1,M+1,2*M+2);
    
    GMatExt = zeros(1,M);
    for a = 1:M
        GMatExt(1,a) = get_G(-0.5,zMat(a),2*pi,L/M); %arbitrary point behind scatterer
    end
    
    N = 1000;
    zExt = linspace(-1,0,N);
    GMatExt2 = zeros(N,M); %Higher resolution grid to account for phase
    for a = 1:N
        for b = 1:M
            GMatExt2(a,b) = get_G(zExt(a),zMat(b),2*pi,L/M); %arbitrary point behind scatterer
        end
    end

    if inputPhase >= 0
        Lphase = inputPhase/(2*pi);
    else
        Lphase = (2*pi+inputPhase)/(2*pi);
    end
    GPhase1 = zeros(1,M);
    maxPoint = -Lphase;
    for j = 1:M
        GPhase1(1,j) = get_G(maxPoint,zMat(j),2*pi,L/M);
    end
    node = maxPoint-0.25;
    GPhase2 = zeros(1,M);
    for j = 1:M
        GPhase2(1,j) = get_G(node,zMat(j),2*pi,L/M);
    end

    GP1Opt = [zeros(M),(1/2)*GPhase1';(1/2)*GPhase1,0];
    GP2Opt = [zeros(M),(1/2)*GPhase2';(1/2)*GPhase2,0];
    
    DMat(:,:,2*M+1) = eye(M);
    DMat(:,:,2*M+2) = 1i*eye(M);
    
    for l = 1:(2*M+2)
        if l <= 2*M
            if(logical(mod(l,2)))
                DMat((l+1)/2,(l+1)/2,l) = 1;
            else
                DMat(l/2,l/2,l) = 1i;
            end
        end
        AInt = DMat(:,:,l)*conj(w)*(GMat-chiInv); 
        %Intermediate AMat element (need to take real part)
        AMat(:,:,l) = (AInt + AInt')/2; %Matrix version of "real(AInt)"
        clear AInt
        BMat(:,:,l) = (1/2)*conj(w)*DMat(:,:,l)*EMat;
        BMatT(:,:,l) = ((1/2)*conj(w)*DMat(:,:,l)*EMat)';
        COpt(:,:,l) = [AMat(:,:,l),BMat(:,:,l);BMatT(:,:,l),0];
    end
    
    phiMax = 1./((GMat(1,:)-chiInv(1,:))*(GMat(1,:)-chiInv(1,:))');
    
    cvx_begin sdp quiet
        variable X(M+1,M+1) hermitian
        minimize(abs(1-trace(GP1Opt*X)));
        subject to
            for a = 1:(2*M+2)
                trace(COpt(:,:,a)*X) == 0;
            end
            X(M+1,M+1) == 1;
            X >= 0;
            trace(GP2Opt*X) == 0;
      cvx_end

    maxP = (1-cvx_optval)^2;

    %Get polarization current from X
    [phiMatInt,~] = eigs(X);
    phiMat = phiMatInt(1:M,1)/phiMatInt(end,1);
    
    Eref = GMatExt2*phiMat;
    [~,maxInd] = max(real(Eref));
    phase = -zExt(maxInd)*2*pi;
    if phase > pi
        phase = phase-2*pi;
    end
    solnMat = X;
end