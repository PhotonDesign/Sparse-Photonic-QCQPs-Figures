function xi_best = Gauss_SDP_dif(X,nRand,Nx,einc,M,rho,chi)
    
    %Make sure X is exactly PSD
    [VX,DX] = eig(X);
    nEigs = length(find(DX > 0));
    X_PSD = zeros(size(X));
    for i = flip((length(X)-nEigs+1):length(X))
        cur_V = VX(:,i);
        cur_eig = DX(i,i);
        X_PSD = X_PSD + cur_eig*(cur_V*cur_V');
    end
    
    %Do 1 run to get initial randomization
    xi_rand = mvnrnd(zeros(length(X_PSD),1),X_PSD)';
    xi_rand = xi_rand(1:(length(xi_rand)-1)) / xi_rand(end);
    xi_rand = xi_rand(1:Nx) + 1j*xi_rand(Nx+1:2*Nx);
    etot = xi_rand + einc;
    errorMat = (M + chi*eye(Nx))*etot + rho;
    error_best = norm(errorMat);
    xi_best = xi_rand;
    clear xi_rand errorMat etot
    
    for i = 2:nRand
        
        xi_rand = mvnrnd(zeros(length(X_PSD),1),X_PSD)';
        xi_rand = xi_rand(1:(length(xi_rand)-1)) / xi_rand(end);
        xi_rand = xi_rand(1:Nx) + 1j*xi_rand(Nx+1:2*Nx);
        etot = xi_rand + einc;
        errorMat = (M + chi*eye(Nx))*etot + rho;

        error = norm(errorMat);

        if error < error_best
            xi_best = xi_rand;
            error_best = error;
        end

        clear xi_rand errorMat error etot
    
    end
    
end