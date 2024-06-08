%% construct Greens function
[X1,Y1] = meshgrid(x,y);
X1_vect = reshape(X1.',[],1);
Y1_vect = reshape(Y1.',[],1);

r2 = zeros(Kx*Ky);

% tic
for ii = 1:Kx*Ky
for jj = 1:Kx*Ky
    
r2(ii,jj) = (X1_vect(ii) - X1_vect(jj))^2 + (Y1_vect(ii) - Y1_vect(jj))^2;

end
end

G0 = (-1i/4).* besselh(0,1,omega*sqrt(r2));
G0(isnan(G0)) = 0; % get rid of singularity
G0_VIE = G0(des_idx_vect,des_idx_vect); % only keep design region for VIE
G0_VIE = (-k^2)*G0_VIE*(h^2);

%% calculate polarization for slab
einc_VIE = einc(Ny+1:My,Nx+1:Mx);  % using same incident field as in diff. calculation (truncated)
einc_vect_VIE = reshape(einc_VIE.',[],1);
Kx_VIE = size(einc_VIE,2);
Ky_VIE = size(einc_VIE,1);

p_slab_vect_VIE = (G0_VIE + (xi*eye(Kx_VIE*Ky_VIE)))\(-einc_vect_VIE);
p_slab_VIE = reshape(p_slab_vect_VIE,[Kx_VIE,Ky_VIE]).';

es_slab_vect_VIE = G0_VIE*p_slab_vect_VIE;
es_slab_VIE = reshape(es_slab_vect_VIE,[Kx_VIE,Ky_VIE]).';

G_tar_VIE = G_tar(des_idx_vect);