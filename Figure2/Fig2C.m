clear
startup

%% user-defined parameter
nD = logspace(1,4, 10);

%% init
TotalDim = nD * 2 + 1;

%% cal
Ny = 13;
maxC_metasurface = 2 * (Ny * 2 + 1) + 1;

maxC_layered = 2 * 3 + 1;

Ny = 5;
maxC_QM = 2 * (Ny * 2 + 1) + 1;

%% result
figure 
set(gcf, 'position', [100, 100, 450, 400])

hold on
plot(nD, TotalDim, 'r-')
plot(nD, maxC_metasurface * ones(size(nD)), 'k--')
plot(nD, maxC_layered * ones(size(nD)), 'k-')
plot(nD, maxC_QM * ones(size(nD)), 'k:')
xlabel('designable degree of freedom')
ylabel('maximal clique size')
set(gca, 'xscale', 'log')
set(gca, 'yscale', 'log')
xticks([200, 400, 800, 1600, 3200])
xlim([169, 3200])
ylim([5e0, 7e3])
yticks(logspace(log10(5), log10(5*2^10), 11))
text(880,1400,'dense matrices', 'color', 'red', 'fontsize', 14)
text(1200,73,'metasurface', 'color', 'black', 'fontsize', 14)
text(1340,32,'waveguide', 'color', 'black', 'fontsize', 14)
text(1080,9.5,'multilayer film', 'color', 'black', 'fontsize', 14)

%% save%% savefig
saveas(gcf, [mfilename, '.eps'], 'epsc')
saveas(gcf, [mfilename, '.png'])

