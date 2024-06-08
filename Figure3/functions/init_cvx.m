function[S, b, c] = init_cvx(Kx,Ky,bg_idx,MM_bg,des_idx,M0,xi,k,einc_vect,lin_bg_const)

if ~exist('lin_bg_const')
    lin_bg_const = 1;
end
    

S = {};
b = {};
c = {};

Mquad = zeros(Kx*Ky);
for ii = bg_idx.'
    
    if lin_bg_const
        
    S(end+1) = {zeros(Kx*Ky)}; % linear constraints for points in the background
    b(end+1) = {MM_bg(ii,:)'};
    c(end+1) = {0};
    
    end
    Mquad = Mquad + MM_bg(ii,:)' * MM_bg(ii,:);
end

S(end+1) = {Mquad}; % summed redundant quadratic constraints
b(end+1) = {zeros(Kx*Ky,1)};
c(end+1) = {0};

for ii = des_idx.'
    Di = zeros(Kx*Ky); % D-matrix for points in the degign region
    Di(ii,ii) = 1;
    S(end+1) = {Di * M0 - (xi'./k^2) * (M0' * Di * M0)};
    b(end+1) = {M0' * Di * einc_vect};
    c(end+1) = {0};
end

% complex to real constraints
[S, b, c] = c2r_con(S, b, c);