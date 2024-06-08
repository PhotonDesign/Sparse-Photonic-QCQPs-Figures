function [maxP_opt,maxP_sim,phiMat,solnMat] = getMaxPsca_rank1(n,L,M,e)
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

    GOpt = [(GMat-GMat')/(2i),zeros(M,1);zeros(M,1)',0];
    
    DMat = zeros(M,M,2*M);
    AMat = zeros(M,M,2*M);
    BMat = zeros(M,1,2*M);
    BMatT = zeros(1,M,2*M);
    COpt = zeros(M+1,M+1,2*M);
    
    for l = 1:(2*M)
        if(logical(mod(l,2)))
            DMat((l+1)/2,(l+1)/2,l) = 1;
        else
            DMat(l/2,l/2,l) = 1i;
        end
        AInt = DMat(:,:,l)*conj(w)*(GMat-chiMat); 
        %Intermediate AMat element (need to take real part)
        AMat(:,:,l) = (AInt + AInt')/2; %Matrix version of "real(AInt)"
        clear AInt
        BMat(:,:,l) = (1/2)*conj(w)*DMat(:,:,l)*EMat;
        BMatT(:,:,l) = ((1/2)*conj(w)*DMat(:,:,l)*EMat)';
        COpt(:,:,l) = [AMat(:,:,l),BMat(:,:,l);BMatT(:,:,l),0];
    end
    
    basMat = eye(M+1);
           
        % Optimization problem, maximize power output
    cvx_begin sdp quiet
        variable X(M+1,M+1) hermitian
        %maximize((conj(w)/2)*trace(GOpt*X)*(L/M) - e*imag(trace(conj(chi)*X)));
        maximize((conj(w)/2)*trace(GOpt*X)*(L/M));
        subject to
            for a = 1:(2*M)
                trace(COpt(:,:,a)*X) == 0;
            end
            X(M+1,M+1) == 1;
            X >= 0;
            %trace(X) <= e;
      cvx_end
    
    maxP_opt = conj(w)*trace(GOpt*X)*(L/M);
    %Get polarization current from X
    [phiMatInt,~] = eigs(X);
    phiMat = phiMatInt(1:M,1)/phiMatInt(end,1);
    maxP_sim = conj(w)*phiMat'*(GMat - GMat')/(2i)*phiMat*L/M;
    %Use 1:M to take off "s"
    solnMat = X;
end