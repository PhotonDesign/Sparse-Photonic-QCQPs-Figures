close all; clc; clear all;
startup

%%
rng(13)
Nm = 10; % number of material patches
L = 8*pi+0.4;
xm = sort(rand(2*Nm,1)) * L;
xm(1) = 0.2;
xm(end) = xm(1)+8*pi;
xm(end-5) = 19.2;
Xm = [xm(1:2:end), xm(2:2:end)]; % coordinates for the material patches
k = 1;
lam = 2*pi/k;

%%
[xx, e_CFEM] = fun_CFEM(802, 1000, 1);
[xx, e_LCFEM] = fun_CFEM(802, 1000, 1/16);
[xx, e_FEM] = fun_FEM(801, 1000);


%%
figure
set(gcf, 'position', [100, 100, 850, 500])

subplot(4, 1, 1)
Xm2 = Xm-xm(1);
for i = 1:Nm
    rectangle('Position',[Xm2(i, 1), 0.001, Xm2(i, 2) - Xm2(i, 1), 0.998], 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', [0.5 0.5 0.5])
end
xlim([0, xm(end)-xm(1)])
ylim([0, 1])

set(gca, 'xtick', [])
set(gca, 'ytick', [])

subplot(4,1,2:4)
hold on
plot(xx / lam, real(e_CFEM) * 2, '-', 'linewidth', 6, 'color', [1 .0 .0])
plot(xx / lam, real(e_LCFEM) * 2, '--', 'linewidth', 4, 'color', [.9 .7 .3])
plot(xx / lam, imag(e_CFEM) * 2, 'linewidth', 6, 'color', [0 0 1])
plot(xx / lam, imag(e_LCFEM) * 2, '--', 'linewidth', 4, 'color', [.3 .7 0.9])
plot(xx / lam, real(e_FEM*2), 'linewidth', 2, 'color', [0 .0 .0]);
plot(xx / lam, imag(e_FEM*2), 'linewidth', 2, 'color', [.0 .0 .0]);
axis([0,4,-0.7,0.7])
yticks([-0.5,0,0.5])
xticks([0 1 2 3 4])
set(gca,'Yticklabel',[]) 
set(gca,'Xticklabel',[]) 






