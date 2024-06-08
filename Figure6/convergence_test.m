close all; clc; clear all;
startup
%% analytical solution
Nx = 15731;

[xx, field] = fun_analytical(Nx);

%% calculate error of CFEM/LCFEM

resolution_list = [102 200 300 400 500 603 700 800 900 1000 2000 3000 4000 5001 6000 6999 8001 9000 9999];
error_list_CFEM = [];
error_list_LCFEM = [];
for resolution = resolution_list
    [~, e_CFEM] = fun_CFEM(resolution, Nx, 1);
    e_CFEM = e_CFEM(1:end-1)*2;
    error_CFEM = sum(abs(e_CFEM-field).^2)*(xx(2)-xx(1));
    [~, e_LCFEM] = fun_CFEM(resolution, Nx, 1/16);
    e_LCFEM = e_LCFEM(1:end-1)*2;
    error_LCFEM = sum(abs(e_LCFEM-field).^2)*(xx(2)-xx(1));
    error_list_CFEM = [error_list_CFEM, error_CFEM];
    error_list_LCFEM = [error_list_LCFEM, error_LCFEM];

end

%% calculate error of FEM


resolution_list2 = [102 202 299 401 500 600 701 801 901 1001 2001 3001 4001 4999 6000 7000 7999 9000 10000];
error_list_FEM = [];
for resolution = resolution_list2
    [~, e_FEM] = fun_FEM(resolution, Nx);
    e_FEM = e_FEM(1:end-1)*2;
    error_FEM = sum(abs(e_FEM-field).^2)*(xx(2)-xx(1));
    error_list_FEM = [error_list_FEM, error_FEM];

end


%%
figure
set(gcf, 'position', [100, 100, 850, 500/2])
loglog(resolution_list2, error_list_FEM,'-o','color',[0 0 0],'LineWidth',3)
hold on
loglog(resolution_list, error_list_CFEM,'-o','LineWidth',3)
loglog(resolution_list, error_list_LCFEM,'--o','LineWidth',3)

axis([100,10000,1e-5,1e1])
yticks([1e-5 1e-3 1e-1 1e1])
