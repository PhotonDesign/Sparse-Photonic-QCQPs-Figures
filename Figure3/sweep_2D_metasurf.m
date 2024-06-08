% lambda = from sweep wrapper
omega = 2*pi/lambda; % assuming c=1;
k = omega; % assuming c=1;

obj = 'focusing';

% chi = from sweep wrapper
xi = -1/chi;

% Ly = from sweep wrapper
% Lx = from sweep wrapper
% h = from sweep wrapper

% %calculation mode
% const_NA_or_d = from sweep wrapper

% d = from sweep wrapper
% NA = from sweep wrapper

% %also do VIE calculation for comparison
% do_VIE = = from sweep wrapper
% 
% % angles to sweep
% Ntheta = from sweep wrapper
% dTheta = from sweep wrapper

plots_on = 0; % plot while running

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Dim = [Lx, Ly]; %[x,y]

Kx = round(Lx/h); %num of x dim grid points
Ky = round(Ly/h); %num of y dim grid points

% thickness = from sweep wrapper
% beta = from sweep wrapper
BC = {{'pml', [thickness,beta]}, {'pml', [thickness,beta]}};
        
Matx = ones(Ky,Kx); %mux, staggered
Maty = ones(Ky,Kx); %muy, staggered
Matz = ones(Ky,Kx); %epsz, not staggered - vacuum background
        
[pmlx, pmly, pmlz] = PML(Dim,h,BC); %get PML
[PML_bool,PML_bool_vect,PMLx_num,PMLy_num] =...
    get_PML_bool(thickness,Kx,Ky,Lx,Ly); %indices of PML regions

x  = linspace(0,Lx-h,Kx);
y  = linspace(0,Ly-h,Ky);

Dim_mat = Dim - 2*(thickness+h);

Lx_mat = Dim_mat(1);
Ly_mat = Dim_mat(2);

x0_mat = 0.5*(Lx-Lx_mat-h); %place dielectric in the center of the simulation
y0_mat = 0.5*(Ly-Ly_mat-h); %place dielectric in the center of the simulation
        
Nx = sum(x <= x0_mat); %last point before the grid enters the design region
Mx = sum(x < x0_mat+Lx_mat); %last point before the grid exits the design region
        
Ny = sum(y <= y0_mat); %last point before the grid enters the design region
My = sum(y < y0_mat+Ly_mat); %last point before the grid exits the design region
        
% get vector indices of designable and non-designable regions
[des_bool,des_bool_vect,bg_idx_vect,des_idx_vect] = get_des_bool(Kx,Ky,Nx,Mx,Ny,My);

% unstructrued dielectric
Matz_slab = Matz;
Matz_slab(Ny+1:My,Nx+1:Mx) = Matz(Ny+1:My,Nx+1:Mx) + chi;

% get differential operator
MM_bg = -Scatt_Maxwell_Operator_Construct...
    (omega,Dim,h,BC,(Matx.*pmlx).',(Maty.*pmly).',(Matz.*pmlz).');
        
MM_slab = -Scatt_Maxwell_Operator_Construct...
    (omega,Dim,h,BC,(Matx.*pmlx).',(Maty.*pmly).',(Matz_slab.*pmlz).');


%% set source (plane wave)      
sx = Nx; % define 0-phase plane

einc_pw_fun = @(xx) repmat(exp(1i*k*(xx - x(sx))),length(y),1); % plane wave
einc = einc_pw_fun(x);
einc_vect = reshape(einc.',[],1);

J_pw_vect = (1i/omega)*MM_bg*einc_vect; %%% need to think about these steps
J_pw = reshape(J_pw_vect,[Kx,Ky]).';

e_slab_vect = MM_slab\(-1i*omega*J_pw_vect); % total field with dielectric slab
e_slab = reshape(e_slab_vect,[Kx,Ky]).';

es_slab_vect = e_slab_vect - einc_vect; % scattered field with dielectric slab

if plots_on
%% plot geometry
figure(1);set(gcf,'Position',[500,300,500,500]);clf;colormap(jet)
subplot 141;imagesc(abs(PML_bool));title('PML');xlim([1,Kx]);ylim([1,Ky])
subplot 142;imagesc(abs(des_bool));title('design');xlim([1,Kx]);ylim([1,Ky])
subplot 143;imagesc(real(J_pw));title('Re(J)');xlim([1,Kx]);ylim([1,Ky])
subplot 144;imagesc(imag(J_pw));title('Im(J)');xlim([1,Kx]);ylim([1,Ky])
plot_set(gcf,12,1);drawnow;
end

if plots_on
%% plot incident field
figure(2);set(gcf,'Position',[1000,300,850,500]);clf;colormap(jet)
subplot 141; imagesc(real(einc));
title('Re(Einc)');xlim([1,Kx]);ylim([1,Ky]);colorbar('southoutside')
subplot 142; imagesc(imag(einc)); 
title('Im(Einc)');xlim([1,Kx]);ylim([1,Ky]);colorbar('southoutside')
subplot 143; imagesc(real(e_slab));
title('Re(E slab)');xlim([1,Kx]);ylim([1,Ky]);colorbar('southoutside')
subplot 144; imagesc(imag(e_slab));
title('Im(E slab)');xlim([1,Kx]);ylim([1,Ky]);colorbar('southoutside')
plot_set(gcf,12,1);drawnow
end

%% set target point
if strcmp(const_NA_or_d,'NA')
d = (Ly_mat/2)/NA*sqrt(1-NA^2);
end

x_tar = x(Mx) + d; 
y_tar = Ly/2;

einc_tar = einc_pw_fun(x_tar);
einc_tar = einc_tar(round(Ky/2)); % incident wave at target point

G_tar = get_G_tar(omega,x,y,x_tar,y_tar);
G_tar(bg_idx_vect) = 0; % only need to integrate in design region
G_tar = (-k^2)*G_tar*(h^2); % Green's function to target point

% full field at target point for a slab using polarization
p_slab = (Matz_slab-1).*e_slab;
p_slab_vect = reshape(p_slab.',[],1);
e_tar = einc_tar + G_tar*p_slab_vect;

theta_slab = angle(e_tar); % use as starting point for angle sweeping

if do_VIE % using VIE
init_focus_VIE;
e_tar_VIE = einc_tar + G_tar_VIE*p_slab_vect_VIE;
end

% %% get objective function and calculate for unstructured waveguide
% Ntheta_tmp = 10;
% theta_list_tmp = linspace(0,2*pi-(2*pi/Ntheta_tmp),Ntheta_tmp); % test theta dependence in target function for slab
% theta_list_tmp = sort(mod(theta_list_tmp - theta_slab,2*pi));
% target_norm = -inf;
% 
% for iii = 1:length(theta_list_tmp)
% theta_tmp = theta_list_tmp(iii);
% 
% [S0,b0,c0] = get_obj_focus(obj,theta_tmp,G_tar,MM_bg,omega,einc_tar);
% fobj_slab(iii) = es_slab_vect'*S0*es_slab_vect + real(b0'*es_slab_vect) + c0;
% 
% if do_VIE % using VIE
% [A_VIE,b_VIE,c_VIE] = get_obj_focus_VIE(obj,theta_tmp,G_tar_VIE,einc_tar);
% fobj_slab_VIE(iii) = p_slab_vect_VIE'*A_VIE*p_slab_vect_VIE...
%     + real(b_VIE'*p_slab_vect_VIE) + c_VIE; 
% end
% 
% end
% 
% if plots_on
% %%
% figure(3);set(gcf,'Position',[1300,600,500,300]);clf;colormap(jet)
% 
% plot(theta_list_tmp/pi,fobj_slab.^2,'ro-');hold on
% plot(theta_list_tmp/pi,abs(e_tar)^2 +0*theta_list_tmp,'k--')
% plot(mod(-theta_slab/pi,2)*[1,1],[0,abs(e_tar)^2],'m:')
% 
% if do_VIE % using VIE
% plot(theta_list_tmp/pi,fobj_slab_VIE.^2,'bs-');hold on
% plot(theta_list_tmp/pi,abs(e_tar_VIE)^2 +0*theta_list_tmp,'k-.')
% lgd = legend('Diff.','Diff.','','VIE','VIE','location','best');
% end
% 
% xlabel('\theta/\pi');ylabel('|E tar|^2');xlim([0,2]);ylim([0,1.5*max(abs(e_tar)^2)]);
% lgd.ItemTokenSize = [12,18];
% 
% title('slab')
% plot_set(gcf,12,1);drawnow
% end


%% calculate bounds
[S,b,c] = init_SparseCoLo(bg_idx_vect,des_idx_vect,MM_bg,xi,k,einc_vect); % get constraints

theta_list = linspace(0,theta_range-(theta_range/Ntheta),Ntheta);
theta_list = -theta_slab + theta_list - mean(theta_list); % calculate bound for different phases (theta)

% init
global Time
fobj_bd = nan(size(theta_list));
time_init_per_theta = nan(size(theta_list));
time_solve_per_theta = nan(size(theta_list));
time_total_per_theta = nan(size(theta_list));

if do_VIE
fobj_bd_VIE = nan(size(theta_list));
time_VIE_solve_per_theta = nan(size(theta_list));
time_VIE_per_theta = nan(size(theta_list));
end

time_Start = tic;

% Ntheta = 0; %only for debugging

for ii = 1:Ntheta
disp(' ');disp('-------------------------------------------------------------------')
disp(['    ','theta: iteration no. ',num2str(ii),' out of ',num2str(Ntheta)]);disp(' ')

time_per_theta_Start = tic;

theta = theta_list(ii);
[S0,b0,c0] = get_obj_focus(obj,theta,G_tar,MM_bg,omega,einc_tar); % get target function
fobj_bd(ii) = complexQCQP_SparseCoLO(S0,b0,c0,S,b,c,'mosek','max');

time_init_per_theta(ii) = Time.SparseCoLO_init;
time_solve_per_theta(ii) = Time.SparseCoLO;
time_total_per_theta(ii) = toc(time_per_theta_Start);

%% compare with a VIE calculation using polarization
if do_VIE   
VIE_tic = tic;

[A_VIE,b_VIE,c_VIE] = get_obj_focus_VIE(obj,theta,G_tar_VIE,einc_tar); % objective function

AA = [A_VIE,1/2*b_VIE;1/2*b_VIE',c_VIE]; % homogenized matrix

BB = @(D) [Re(D*omega'*(G0_VIE+(xi*eye(Kx_VIE*Ky_VIE)))),1/2*omega'*D*einc_vect_VIE ;...
    1/2*omega*einc_vect_VIE'*D' , 0]; % constraints

D = get_all_D(Kx_VIE*Ky_VIE); % all D matrices

VIE_solve_tic = tic;

[fobj_bd_VIE(ii),~] = get_bound_VIE(D,AA,BB,1); % last argument sets cvx to quiet

time_VIE_solve_per_theta(ii) = toc(VIE_solve_tic);
time_VIE_per_theta(ii) = toc(VIE_tic);

end

end

time_end = toc(time_Start);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%