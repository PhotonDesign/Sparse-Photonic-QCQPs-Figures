close all; clc; clear all;
set(0,'DefaultLineLineWidth', 1);
set(0,'DefaultAxesTitleFontWeight','normal');

% cvx_solver mosek
% cvx_save_prefs
%addpath(genpath('C:\Program Files\Mosek'))
addpath(genpath('functions'))
%addpath(genpath('../'))

unit = 1e6; % scale length

eps0 = 8.854e-12/unit;
mu0 = pi*4e-7/unit;
c_light = 1/sqrt(eps0*mu0);

eps_r = 1; % background material
lambda = 1.55;
omega = 2*pi*c_light/lambda;
k = omega/c_light;

Dim = [2, 0.03];
h = 0.0075; % discretization step size

%design region
n_mat = sqrt(12);
L_mat = 0.7;

chi = (n_mat^2)-1; % chi of dielectric
xi = -1/chi;

obj = 'scattering';
% obj = 'extinction';
obj = 'integrate full field';
% obj = 'full field at center';
% obj = 'full field at center line';
% obj = 'integrate full field rect';

L_mat_LIST = (0.02:0.04:.9); % to plot geometry no need to sweep 

if length(L_mat_LIST)==1
    plots_on = 1; % plot while running
else
    plots_on = 0;
end

%%
time_tot_start = tic;
for mm = 1:length(L_mat_LIST)
disp(' ');disp('----------------------');
disp(['L_mat: iteration no. ',num2str(mm),' out of ',num2str(length(L_mat_LIST))])
disp(' ');
L_mat = L_mat_LIST(mm);

%%
Lx = Dim(1);
Ly = Dim(2);

Kx = round(Lx/h); %num of x dim grid points
Ky = round(Ly/h); %num of y dim grid points

% x  = linspace(0,Lx-h,Kx);
% y  = linspace(0,Ly-h,Ky);
x  = linspace(0,Lx,Kx);
y  = linspace(0,Ly,Ky);

thickness = 0.5; %pml thickness
beta = 3e7/unit; %pml strength
BC = {{'pml', [thickness,beta]}, {'periodic', [thickness,beta]}};

Matx = ones(Ky,Kx); %mux, staggered
Maty = ones(Ky,Kx); %muy, staggered
Matz = eps_r*ones(Ky,Kx); %epsz, not staggered

% get PML
[pmlx, pmly, pmlz] = PML(Dim,h,BC);
%indices of PML regions
PMLx_num = round(thickness*Kx/Lx);
PML_bool = zeros(Ky,Kx);
PML_bool(:,1:PMLx_num) = 1;
PML_bool(:,end-PMLx_num+1:end) = 1;
PML_bool_vect = reshape(PML_bool.',[],1);

%% design region
Dim_mat = [L_mat,Ly]; % size of design region (x,y)

Lx_mat = Dim_mat(1);
Ly_mat = Dim_mat(2);

% x0_mat = 0.5*(Lx-Lx_mat-h); %place dielectric in the center of the simulation
% y0_mat = 0.5*(Ly-Ly_mat-h); %place dielectric in the center of the simulation

x0_mat = 0.5*(Lx-Lx_mat); %place dielectric in the center of the simulation
y0_mat = 0.5*(Ly-Ly_mat); %place dielectric in the center of the simulation

Nx = sum(x <= x0_mat); %last point before the grid enters the design region
Mx = sum(x < x0_mat+Lx_mat); %last point before the grid exits the design region

Ny = 0; sum(y <= y0_mat); %last point before the grid enters the design region
My = Ky; sum(y < y0_mat+Ly_mat); %last point before the grid exits the design region

% get vector indices of designable and non-designable regions
[des_bool,des_bool_vect,bg_idx_vect,des_idx_vect] = get_des_bool(Kx,Ky,Nx,Mx,Ny,My);

if isempty(des_idx_vect)
    disp('Design region is empty!');
    return;
end

% source
sx = round(0.6*Kx/(Lx-h));
Jsource = -1i*omega*mu0/(Ly*h);

Source = zeros(Ky,Kx);
Source(:,sx) = Jsource;
Source_vect = reshape(Source.',[],1);

if strcmp(obj,'full field at center')% target point at center of material
Target_point = zeros(Ky,Kx);
Target_point(round(Ky/2),round(Kx/2)) = 1;
Target = Target_point;

elseif strcmp(obj,'full field at center line') % target line at center of material
Target_line = zeros(Ky,Kx);
Target_line(:,round(Kx/2)) = 1;
Target = Target_line;

elseif strcmp(obj,'integrate full field rect') % target rectangle
rect_L = 0.25;
rect_center = Lx/2;
    
Target_rect = zeros(Ky,Kx);
Target_rect(:,round((rect_center-rect_L/2)/h):round((rect_center+rect_L/2)/h)) = 1;
Target = Target_rect;

else
Target = 0;
end

Target_vect = reshape(Target.',[],1);
target_idx = find(Target_vect);

if plots_on
%% plot geometry
figure(1);set(gcf,'Position',[1000,300,500,600]);clf;
subplot 411; imagesc(PML_bool);axis equal;xlim(inf*[-1,1]);ylim(inf*[-1,1]);title('PML')
subplot 412; imagesc(abs(Source));axis equal;xlim(inf*[-1,1]);ylim(inf*[-1,1]);title('|source|')
subplot 413; imagesc(des_bool);axis equal;xlim(inf*[-1,1]);ylim(inf*[-1,1]);title('design region')
subplot 414; imagesc(Target);axis equal;xlim(inf*[-1,1]);ylim(inf*[-1,1]);title('target')
plot_set(gcf,12,1);colormap(flip(pink));drawnow;
end

%% get operators
MM_bg = -Scatt_Maxwell_Operator_Construct...
    (omega/c_light,Dim,h,BC,(Matx.*pmlx).',(Maty.*pmly).',(Matz.*pmlz).');

% unstructrued dielectric
Matz_slab = Matz;
Matz_slab(Ny+1:My,Nx+1:Mx) = Matz(Ny+1:My,Nx+1:Mx) + chi;

MM_slab = -Scatt_Maxwell_Operator_Construct...
    (omega/c_light,Dim,h,BC,(Matx.*pmlx).',(Maty.*pmly).',(Matz_slab.*pmlz).');

% incident field
einc_vect = MM_bg\Source_vect;
einc = reshape(einc_vect,[Kx,Ky]).';

% dielectric slab field
e_slab_vect = MM_slab\Source_vect;
e_slab = reshape(e_slab_vect,[Kx,Ky]).';

% dielectric slab  scattered field
es_slab_vect = e_slab_vect - einc_vect;
es_slab = e_slab - einc;

if plots_on
%% plot incident and slab fields
figure(2);set(gcf,'Position',[1000,500,600,400]);clf;

subplot 311;imagesc(real(einc));hold on
axis equal;plot([sx,sx],[1,Ky],'wo--')
xlim(inf*[-1,1]);ylim(inf*[-1,1]);title('Re(Einc)');colorbar('southoutside')

subplot 312; imagesc(real(e_slab));hold on
axis equal;plot([Nx+1,Nx+1],[1,Ky],'kx--',[Mx,Mx],[1,Ky],'kx--');
xlim(inf*[-1,1]);ylim(inf*[-1,1]);title('Re(E slab)');colorbar('southoutside')

subplot 313; imagesc(real(es_slab));
axis equal;xlim(inf*[-1,1]);ylim(inf*[-1,1]);title('Re(Es slab)');colorbar('southoutside')

plot_set(gcf,12,1);colormap(jet);drawnow;
end

%% get target function and calculate for slab field

% Wm = get_trapz_mat(Kx,Ky,Nx,Mx,Ny,My)*(h^2); % 2D trapezoidal integration matrix for design region
Wm = diag(des_bool_vect)*(h^2); % multiply matrices without quadrature correction

if strcmp(obj,'scattering')
[S0,b0,c0] = get_obj_scat(obj,Wm,MM_bg,omega,k,eps0,einc_vect); % get target function    
end

if strcmp(obj,'extinction')
[S0,b0,c0] = get_obj_ext(obj,Wm,MM_bg,omega,k,eps0,einc_vect); % get target function    
end

if strcmp(obj,'integrate full field')
[S0,b0,c0] = get_obj_int_e(obj,Wm,einc_vect); % get target function
end

if strcmp(obj,'full field at center')  
[S0,b0,c0] = get_obj_e_cent(obj,target_idx,einc_vect); % get target function
end

if strcmp(obj,'integrate full field rect')
Wm_rect = zeros(Kx*Ky,1);
Wm_rect(target_idx) = h^2; % no trapezoidal rule
Wm_rect = diag(Wm_rect);

[S0,b0,c0] = get_obj_int_e(obj,Wm_rect,einc_vect); 
end

if strcmp(obj,'full field at center line')
% Wm_line = zeros(Kx*Ky,1);
% Wm_line(target_idx(1)) = (1/2)*h;
% Wm_line(target_idx(end)) = (1/2)*h;
% Wm_line(target_idx(2:end-1)) = h;

Wm_line = zeros(Kx*Ky,1); % no trapezoidal rule
Wm_line(target_idx) = h; % no trapezoidal rule
Wm_line = diag(Wm_line);

[S0,b0,c0] = get_obj_int_e(obj,Wm_line,einc_vect); % get target function
end

f_slab = es_slab_vect'*S0*es_slab_vect + real(b0'*es_slab_vect) + c0;

%% calculate bound
time_start = tic;

[S,b,c] = init_SparseCoLo(bg_idx_vect,des_idx_vect,MM_bg,xi,k,einc_vect); % get constraints
[f_bd] = complexQCQP_SparseCoLO(S0,b0,c0,S,b,c,'mosek','max');

time_end = toc(time_start);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

%% COMSOL
k_COMSOL = k; % need to change units??

mesh_res = 2; % 1=extremely fine -- 9=extremely coarse

% get model for incident field
model_air = quasi_1D(1,L_mat,mesh_res);
mphsave(model_air,'Matlab generated Einc'); % saves simulation .mph
s = 'sol3';
xminfo = mphxmeshinfo(model_air,'soltag',s);
coords = xminfo.dofs.coords;

q_air = mphmatrix(model_air,s,'Out',{'L','K'});

if plots_on
%%
figure(11); % plot directly from COMSOL
subplot 221
mphplot(model_air,'pg3','rangenum',1)
axis tight;
subplot 222
mphplot(model_air,'pg2','rangenum',1)
axis tight;
end

%% get model for slab field
model_glass = quasi_1D(n_mat,L_mat,mesh_res);
mphsave(model_glass,'Matlab generated Eslab'); % saves simulation .mph
q_slab = mphmatrix(model_glass,s,'Out',{'L','K'});

if plots_on
%%
subplot 223  % plot directly from COMSOL
mphplot(model_glass,'pg3','rangenum',1)
axis tight;
subplot 224
mphplot(model_glass,'pg2','rangenum',1)
axis tight;
end

if plots_on
%%
figure(12); % plot mesh directly from COMSOL
mphmesh(model_glass);

daspect([1 .1 1])
axis tight;
end

%% organize matrices and sources
COMSOL_M_bg = conj(q_air.K);
COMSOL_M_slab = conj(q_slab.K);
COMSOL_J = sparse(conj(q_air.L));

COMSOL_unit = 1e-6;
COMSOL_J = COMSOL_J*COMSOL_unit;

xx = coords(1,:).';
yy = coords(2,:).';
sz = 15; % for plotting

if plots_on
%%
figure(13);
subplot 311;
scatter(xx,yy,sz,abs(COMSOL_J),'filled')
axis equal;axis tight;colorbar;title('|J|');
end

% find indices of design region and background region
COMSOL_des_idx = find(xx>=-(L_mat*COMSOL_unit)/2 & xx<=(L_mat*COMSOL_unit)/2);
COMSOL_bg_idx = (setdiff(1:length(xx),COMSOL_des_idx)).';

[~,target_idx_C] = min(abs(xx)+abs(yy)); % find point at center (0,0)

if strcmp(obj,'integrate full field rect')
target_idx_C = find(xx>=-(rect_L*COMSOL_unit)/2 & xx<=(rect_L*COMSOL_unit)/2);
end

if strcmp(obj,'full field at center line')
[~,target_idx_C] = min(abs(xx)+abs(yy)); % find point at center (0,0)
target_idx_C = find(xx==xx(target_idx_C)); % define vertical line
end


if plots_on
%%
subplot 312; % plot design region
des_vect = 0*COMSOL_J;
des_vect(COMSOL_des_idx) = 1;
scatter(xx,yy,sz,des_vect,'filled')
axis equal;axis tight;colorbar
title('design region')

subplot 313; % plot target point
tar_vect = 0*COMSOL_J;
tar_vect(target_idx_C) = 1;
scatter(xx,yy,sz,tar_vect,'filled')
axis equal;axis tight;colorbar
title('target')
end



if ~true
    %%
figure(68146)   

ssz = 5;

scatter(xx,yy,ssz,des_vect,'filled')
axis equal;axis tight;colorbar
title('design region')

colormap(jet)

daspect([1 .1 1])
    
end


%% calculate fields from matrices
COMSOL_einc = COMSOL_M_bg\COMSOL_J;
COMSOL_e_slab = COMSOL_M_slab\COMSOL_J;
COMSOL_es_slab = COMSOL_e_slab - COMSOL_einc;

if plots_on
%%
figure(15);
subplot 311; scatter(xx,yy,sz,real(COMSOL_einc),'filled')
axis equal;axis tight;colorbar('southoutside');title('Re(Einc)');
subplot 312; scatter(xx,yy,sz,real(COMSOL_e_slab),'filled')
axis equal;axis tight;colorbar('southoutside');title('Re(E slab)');
subplot 313; scatter(xx,yy,sz,real(COMSOL_es_slab),'filled')
axis equal;axis tight;colorbar('southoutside');title('Re(Es slab)');
colormap(jet)
end

%% interpolate to get smoother plots and field values
F_inc = scatteredInterpolant(xx,yy,full(COMSOL_einc));
F_slab = scatteredInterpolant(xx,yy,full(COMSOL_e_slab));
F_es_slab = scatteredInterpolant(xx,yy,full(COMSOL_es_slab));

interp_res = 1e-9; % resolution for interpolation

x_hi_res = (min(xx)+interp_res):interp_res:(max(xx)-interp_res);
y_hi_res = (min(yy)+interp_res):interp_res:(max(yy)-interp_res);
[XF,YF] = meshgrid(x_hi_res,y_hi_res);

V_inc = F_inc(XF,YF);
V_slab = F_slab(XF,YF);
V_es_slab = F_es_slab(XF,YF);

if plots_on
%%
figure(151);
subplot 311; scatter(XF(:),YF(:),sz,real(V_inc(:)),'filled');
axis equal;axis tight;colorbar('southoutside');title('Re(Einc)');
subplot 312; scatter(XF(:),YF(:),sz,real(V_slab(:)),'filled');
axis equal;axis tight;colorbar('southoutside');title('Re(E slab)');
subplot 313; scatter(XF(:),YF(:),sz,real(V_es_slab(:)),'filled');
axis equal;axis tight;colorbar('southoutside');title('Re(Es slab)');
colormap(jet)
end

%%
if plots_on
%% look at COMSOL matrices
figure(16);
th = 1e-8;
subplot 321;spy(COMSOL_M_bg);title('M bg')
subplot 322;spy(COMSOL_M_slab);title('M slab')
subplot 323;spy(abs(COMSOL_M_bg - COMSOL_M_slab));title('|M bg - M slab|')
subplot 324;spy(abs(COMSOL_M_bg - COMSOL_M_slab)>th)
nnzrows = length(find(sum(abs(COMSOL_M_bg - COMSOL_M_slab)>th)));
title(['th = ',num2str(th),', nnz rows = ',num2str(nnzrows)])

thlist = logspace(-10,-7,50);
for iii = 1:length(thlist)
ppp(iii) = nnz(abs(COMSOL_M_bg - COMSOL_M_slab)>thlist(iii));
qqq(iii) = length(find(sum(abs(COMSOL_M_bg - COMSOL_M_slab)>thlist(iii))));
end

subplot 325;semilogx(thlist,ppp,'x-');xlabel('threshold');ylabel('nonzeros')
subplot 326;semilogx(thlist,qqq,'x-');xlabel('threshold');ylabel('nonzero rows');hold on
semilogx(thlist,0*thlist+length(COMSOL_des_idx),'-.')
drawnow
end

%% get target function and calculate for slab field

Wm_C = get_COMSOL_tri_trapz_mat(xminfo,COMSOL_des_idx);
Wm_C = Wm_C./(COMSOL_unit^2);

if strcmp(obj,'scattering')
[S0_C,b0_C,c0_C] = get_obj_scat(obj,Wm_C,COMSOL_M_bg,omega,k_COMSOL,eps0,COMSOL_einc); % get target function
end

if strcmp(obj,'extinction')
[S0_C,b0_C,c0_C] = get_obj_ext(obj,Wm_C,COMSOL_M_bg,omega,k_COMSOL,eps0,COMSOL_einc); % get target function
end

if strcmp(obj,'integrate full field')
[S0_C,b0_C,c0_C] = get_obj_int_e(obj,Wm_C,COMSOL_einc); % get target function
end

if strcmp(obj,'full field at center')  
[S0_C,b0_C,c0_C] = get_obj_e_cent(obj,target_idx_C,COMSOL_einc); % get target function
end

if strcmp(obj,'integrate full field rect')
Wm_rect_C = get_COMSOL_tri_trapz_mat(xminfo,target_idx_C);
Wm_rect_C = Wm_rect_C./(COMSOL_unit^2);    

[S0_C,b0_C,c0_C] = get_obj_int_e(obj,Wm_rect_C,COMSOL_einc); % get target function   
end


if strcmp(obj,'full field at center line')
    
tmp = yy(target_idx_C);
[t,tt] = sort(tmp);
d = [diff(t);0] +[0;diff(t)];
     
Wm_line_C = zeros(size(COMSOL_einc,1),1);
Wm_line_C(target_idx_C(tt)) = d./COMSOL_unit;
Wm_line_C = diag(Wm_line_C);

[S0_C,b0_C,c0_C] = get_obj_int_e(obj,Wm_line_C,COMSOL_einc); % get target function
end

f_slab_C = COMSOL_es_slab'*S0_C*COMSOL_es_slab + real(b0_C'*COMSOL_es_slab) + c0_C;

%% calculate bounds using COMSOL matrices
time_start_C = tic;

[S_C,b_C,c_C] = init_SparseCoLo(COMSOL_bg_idx,COMSOL_des_idx,COMSOL_M_bg,xi,k_COMSOL,COMSOL_einc); % get constraints
[f_bd_C] = complexQCQP_SparseCoLO(S0_C,b0_C,c0_C,S_C,b_C,c_C,'mosek','max');

time_end_C = toc(time_start_C);

%%
time_end_list(mm) = time_end;
f_slab_list(mm) = f_slab;
f_bd_list(mm) = f_bd;
time_end_C_list(mm) = time_end_C;
f_slab_C_list(mm) = f_slab_C;
f_bd_C_list(mm) = f_bd_C;
Npoints(mm) = Kx*Ky;
Npoints_C(mm) = length(xx);
end

time_tot = toc(time_tot_start);

disp('DONE!')
disp(['total time: ',num2str(time_tot)])

if length(L_mat_LIST)==1
%% print results
disp('-----------------------------------------')
disp(obj);disp(' ')

disp('FDFD')
disp(['    time = ',num2str(time_end),' s']);
disp(['    slab = ',num2str(f_slab,'%.4g')]);
disp(['    bound = ',num2str(f_bd,'%.4g')]);

disp('COMSOL')
disp(['    time = ',num2str(time_end_C),' s']);
disp(['    slab = ',num2str(f_slab_C,'%.4g')]);
disp(['    bound = ',num2str(f_bd_C,'%.4g')]);
end


if length(L_mat_LIST)>1
%% plot sweep
figure(10000);set(gcf,'Position',[900,300,950,350]);clf;

subplot 131
plot(L_mat_LIST,time_end_list,'x-',L_mat_LIST,time_end_C_list,'s-');
xlabel('L des');ylabel('time (s)');axis square
legend('FDFD','COMSOL','location','best');

subplot 132
plot(L_mat_LIST,real(f_slab_list),'x-',L_mat_LIST,real(f_slab_C_list),'s-');
hold on
plot(L_mat_LIST,real(f_bd_list),'x--',L_mat_LIST,real(f_bd_C_list),'s--');
xlabel('L des');ylabel(obj);axis square
legend('FDFD slab','COMSOL slab','FDFD bound','COMSOL bound','location','best');

subplot 133
plot(L_mat_LIST,Npoints,'x-',L_mat_LIST,Npoints_C,'s-');
xlabel('L des');ylabel('no. of points in calculation');axis square
legend('FDFD','COMSOL','location','best');

plot_set(gcf,12,1);
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

return

%% debugging

%% check constraints
es_slab_vect'*S0*es_slab_vect+real(b0'*es_slab_vect)+c0

Ncon = length(S);
for jj = 1:Ncon
    err(jj) = es_slab_vect'*S{jj}*es_slab_vect+real(b{jj}'*es_slab_vect)+c{jj};
end

figure(4);set(gcf,'Position',[1100,600,500,300]);clf;colormap(jet)
subplot 221; plot(1:Ncon,abs(err),'.')
xlim([1,Ncon]);xlabel('real constraint no.');ylabel('|error|')

th = 1e-12;
errfull = err(1:2:end) + err(2:2:end);
errfull_idx = find(abs(errfull)>th);

subplot 222; plot(1:Ncon/2,abs(errfull),'.'); hold on
plot(errfull_idx,abs(errfull(errfull_idx)),'.');
xlabel('complex constraint no.');ylabel('|error|')
title(['th = ',num2str(th)]);

allidx = [bg_idx_vect;des_idx_vect];
allidx_x = 1:length(allidx);
subplot 223; plot(allidx,'.'); hold on
plot(allidx_x(errfull_idx(1:end-1)),allidx(errfull_idx(1:end-1)),'.'); 
xlabel('complex constraint no.');ylabel('spatial index')

debug_vect = zeros(Kx*Ky,1);
debug_vect(allidx(errfull_idx(1:end-1))) = 1;
debug = reshape(debug_vect,[Kx,Ky]).';
subplot 224; imagesc(debug); colormap(gray)

plot_set(gcf,12,1);

%% debug: analyze polarization currents
%%FDFD
p_slab1 = (Matz_slab-1).*e_slab*eps0;
p_slab1_vect = reshape(p_slab1.',[],1);

p_slab2_vect = (-eps0/k^2)*MM_bg*es_slab_vect;
p_slab2 = reshape(p_slab2_vect,[Kx,Ky]).';

figure
plot(real(p_slab1_vect),'x-'); hold on
plot(real(p_slab2_vect),'s-');

figure 
imagesc(abs(p_slab1));

figure 
imagesc(abs(p_slab2));

%%COMSOL
chi_C = 0*xx;
chi_C(COMSOL_des_idx) = chi;

figure 
scatter(xx,yy,sz,chi_C,'filled');

p_slab1_C = chi_C.*COMSOL_e_slab*eps0;
p_slab2_C = (-eps0/k_COMSOL^2)*COMSOL_M_bg*COMSOL_es_slab;

figure
plot(real(p_slab1_C),'x-'); hold on
figure
plot(real(p_slab2_C),'s-');

figure 
scatter(xx,yy,sz,abs(p_slab1_C),'filled');

figure 
scatter(xx,yy,sz,abs(p_slab2_C),'filled');
