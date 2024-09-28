function[fval,dfdx,fcval,dfcdx,r] = sim_dir_adj(alpha,x,Nx,dx,designInd,n,k,rTar,beta)
chi = chi_from_alpha(alpha,n,beta,Nx,designInd);
[etot,~,einc,r,~] = fdfd_1d(x,dx,chi,k);
fval = abs(r-rTar)^2;
pAdj = 2*conj(r-rTar) / einc(3);
[eadj] = fdfd_1d(x,dx,chi,k,pAdj);

chi0 = n^2 - 1;
dfdx = real((chi0+1i*beta*(1-2*alpha)) .* etot(designInd) .* eadj(designInd));
dfdx = dfdx(:);

% one (trivial) constraint for mma code
fcval=-40+1e-3*sum(alpha);
dfcdx=1e-3*ones(1,length(designInd));
dfcdx=dfcdx(:).';
end