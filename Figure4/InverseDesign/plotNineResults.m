
M = cell(9,1);
xmax = 0;
fcmax = 0;
for i=1:9
    M{i} = load(['data',num2str(i),'.mat']);
    xmax = max([xmax, max(sum(M{i}.xHist,2))]);
    fcmax = max([fcmax, max(max(M{i}.fc))]);
end

figure(1);
xv = linspace(0,1.02*xmax,300);
for i=1:9
    epsr = plotDesign(M{i}.xHist(end,:),xv,3,0);
    sbp = subplot(9,1,i);
    imagesc(epsr)
    colormap gray
    axis off
end

figure(2);
xlim([1 fcmax]);
ylim([0 1]);
yticks(0:0.1:1);
yticklabels({'0','','','','','0.5','','','','','1'})
xlabel('# Simulations');
ylabel('Reflectivity (%)');
for i=1:9
    plot(M{i}.fc,100*M{i}.mHist); hold on;
end
legend({'1','2','3','4','5','6','7','8','9'})