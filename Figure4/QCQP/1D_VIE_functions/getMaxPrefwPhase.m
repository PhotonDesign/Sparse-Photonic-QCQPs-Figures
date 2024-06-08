function [maxP,phiMat,Eref,phase,DInds,solnMat] = getMaxPrefwPhase(n,L,M,inputPhase)
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
    chiMat = (1/chi)*eye(M);
    GMat = zeros(M,M);
    for a = 1:M
        for b = 1:M
            GMat(a,b) = get_G(zMat(a),zMat(b),2*pi,L/M);
        end
    end
   
    DMat(:,:,1) = eye(M);
    DMat(:,:,2) = 1i*eye(M);
    
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
    
    for l = 1:M
        if l > 2
            DMat(:,:,l) = w*diag(diag(phiMatOpt(:,:,l-2)*phiMatOpt(:,:,l-2)'*(GMat-chiMat)'+phiMatOpt(:,:,l-2)*EMat'));
            %Gets diagonal elements of matrix in () puts them in
            %diagonal matrix
        end
        AInt = DMat(:,:,l)*conj(w)*(GMat-chiMat); 
        %Intermediate AMat element (need to take real part)
        AMat(:,:,l) = (AInt + AInt')/2; %Matrix version of "real(AInt)"
        clear AInt
        BMat(:,:,l) = (1/2)*conj(w)*DMat(:,:,l)*EMat;
        BMatT(:,:,l) = ((1/2)*conj(w)*DMat(:,:,l)*EMat)';
        COpt(:,:,l) = [AMat(:,:,l),BMat(:,:,l);BMatT(:,:,l),0];
        if l > 1
            % Optimization problem, maximize power output
            cvx_begin sdp quiet
                variable X(M+1,M+1) hermitian
                minimize(abs(1-trace(GP1Opt*X)));
                subject to
                    for a = 1:l
                        trace(COpt(:,:,a)*X) == 0;
                    end
                    X(M+1,M+1) == 1;
                    X >= 0;
                    trace(GP2Opt*X) == 0;
              cvx_end
                
            maxPs(l-1) = (1/2)*(1-cvx_optval)^2;
            
            %Get polarization current from X
            [phiMatInt,~] = eigs(X);
            phiMatOpt(:,:,l-1) = phiMatInt(1:M,1)/phiMatInt(end,1);
            %Use 1:M to take off "s"
        end

        if l > 2 && (abs(maxPs(l-1) - maxPs(l-2))/maxPs(l-1) < 0.01)
            break
        end
    end
    maxP = 2*maxPs(l-1);
    phiMat = phiMatOpt(:,:,l-1);
    Eref = real(GMatExt2*phiMat);
    [~,maxInd] = max(Eref);
    phase = -zExt(maxInd)*2*pi;
    if phase > pi
        phase = phase-2*pi;
    end
    DInds = l;
    solnMat = X;
end