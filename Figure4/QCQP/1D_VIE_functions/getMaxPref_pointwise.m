function [maxP,phiMat,X] = getMaxPref_pointwise(n,L,M)
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
    
    GMatExt = zeros(1,M);
    for a = 1:M
        GMatExt(1,a) = get_G(-0.5,zMat(a),2*pi,L/M); %arbitrary point behind scatterer
    end
    

    GOpt = [GMatExt'*GMatExt,zeros(M,1);zeros(M,1)',0];
    
    DMat = zeros(M,M,2*M+2);
    AMat = zeros(M,M,2*M+2);
    BMat = zeros(M,1,2*M+2);
    BMatT = zeros(1,M,2*M+2);
    COpt = zeros(M+1,M+1,2*M+2);
    
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
    
            % Optimization problem, maximize power output
    cvx_begin sdp quiet
        variable X(M+1,M+1) hermitian
        maximize((1/2)*trace(GOpt*X));
        subject to
            for a = 1:(2*M+2)
                trace(COpt(:,:,a)*X) == 0;
            end
        X(M+1,M+1) == 1;
        X >= 0;     
    cvx_end

    maxP = 2*cvx_optval;

    %Get polarization current from X
    [phiMatInt,~] = eigs(X);
    phiMat = phiMatInt(1:M,1)/phiMatInt(end,1);
            %Use 1:M to take off "s"
    
    X = X;
end