startup

N = 40;
M = cell(N,1);
xmax = 0;
fcmax = 0;
mValsS = zeros(N,1);
mVals = zeros(N,1);
fcVals = zeros(N,1);
for i=1:N
    M{i} = load(['HundredResultsInteriorPoint/data',num2str(i),'.mat']);
%     M{i} = load(['HundredResultsGradDescent/data',num2str(i),'.mat']);
    xmax = max([xmax, max(sum(M{i}.xHist,2))]);
    fcmax = max([fcmax, max(max(M{i}.fc))]);
    mValsS(i) = 100*min(M{i}.mHist);
    mVals(i) = 100*M{i}.mHist(end);
    fcVals(i) = M{i}.fc(end);
end

% % Rbnds = [87.8 87.9; 76.6 76.7; ]
% % 
% figure;
% xv = linspace(0,1.02*xmax,300);
% rng(11);
% ind = zeros(10,1);
% ind(1:8) = round(N*rand(8,1));
% j = find(mVals>85);
% ind(9) = j(1);
% j = find((mVals<80).*(mVals>75));
% ind(10) = j(1);
% for i=1:10
%     ii = ind(i);
%     epsr = plotDesign(M{ii}.xHist(end,:),xv,3,0);
%     sbp = subplot(10,1,i);
%     imagesc(epsr)
%     colormap gray
%     axis off
% end

%% figures
figure;
for i=1:N
     plot(M{i}.fc,100*M{i}.mHist); hold on;
%    plot(M{i}.fc(2:end),100*M{i}.mHist); hold on;
end
% legend({'1','2','3','4','5','6','7','8','9'})
xlim([0,200])

ylim([0 100]);
yticks(0:20:100);
xlabel('# Simulations');
ylabel('Reflectivity (%)');
%%
figure; 
histogram(mValsS,'FaceColor',0.8*[1 1 1]); 
hold on; 
histogram(mVals);
xlim([0 100])
xlabel('Reflectivity (%)')
ylabel('Optimal designs')
hold on; plot(97.5,1,'b*','MarkerSize',10)