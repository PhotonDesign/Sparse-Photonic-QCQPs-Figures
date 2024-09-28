
function[chi] = chi_from_alpha(alpha,n,beta,Nx,designInd)
chi0 = n^2 - 1;
chi = zeros(Nx,1);
chi(designInd) = chi0*alpha + 1i * beta * alpha.*(1-alpha);
end