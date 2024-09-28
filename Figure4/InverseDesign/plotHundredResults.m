
N = 100;
folder = 'MMAHundredResults';

% N = 50;
% folder = 'MMAFiftyResultsNearZeroStart';
% folder = 'MMAFiftyResultsRandStart';

k = 2*pi;
xDesignStart = 0;
L = 4; % size in wavelengths
Nx = 100;
[x,dx,Nx,designInd,Ndesign] = set_fdfd_grid(L,Nx,xDesignStart);
xDesign = x(designInd);

M = cell(N,1);
xmax = 0;
fcmax = 0;
mValsS = zeros(N,1);
mVals = zeros(N,1);
fcVals = zeros(N,1);
for i=1:N
    M{i} = load([folder,'/data',num2str(i),'.mat']);
    xmax = max([xmax, max(sum(M{i}.xHist,2))]);
    fcmax = max([fcmax, max(max(M{i}.fc))]);
    mValsS(i) = 100*min(M{i}.mHist);
    mVals(i) = 100*M{i}.mHist(end);
    fcVals(i) = M{i}.fc(end);
end

% % % Figure 4(b) left side
figure;
rng(15); % generate same random numbers (to reproduce figure)
ind = zeros(10,1);
ind(1:8) = round(N*rand(8,1));
j = find(mVals>97.5);
ind(9) = j(1);
% j = find((mVals<80).*(mVals>75));
ind(10) = j(2);
for i=1:10
    ii = ind(i);
    alpha = M{ii}.xHist(end,:); 
    alpha(alpha<0.5)=0;
    alpha(alpha>0.5)=1;
    alpha = repmat(alpha,3,1);
    sbp = subplot(10,1,i);
    imagesc(alpha)
    hAxes = gca;
    colormap(hAxes, [1 1 1; 0 0 0]);
    axis off
end
mVals(ind) % prints the efficiency values
% 
% % % Compare the spacings in an initial topology-optimization step versus
% % % an optimal design
% figure;
% rng(15); % generate same random numbers (to reproduce figure)
% ind = zeros(10,1);
% ind(1:8) = round(N*rand(8,1));
% ind = ind(1);
% iter = 2;
% alpha = M{ii}.xHist(iter,:);
% sbp = subplot(2,1,1);
% imagesc(xDesign,[-1 0 1], alpha)
% hAxes = gca;
% colormap(hAxes, [1 1 1; 0 0 0]);
% set(hAxes,'ytick',[]);
% d = load('../QCQP-Result/Rank1_MM_4lambda.mat');
% alpha = real(d.chi_dif_b(:).' / max(d.chi_dif_b(:))); % convert chi values back to density
% alpha = repmat(alpha,3,1);
% sbp = subplot(2,1,2);
% imagesc(d.x/2/pi, [-1 0 1], alpha)
% hAxes = gca;
% colormap( hAxes , [1 1 1; 0 0 1] )
% set(hAxes,'ytick',[]);

% % % Figure 4(b) right-hand side
% figure;
% for i=1:N
%     plot(M{i}.fc,100*M{i}.mHist,'color',[0.7 0.7 0.7]); hold on;
%     plot(M{i}.fc,100*M{i}.mHist); hold on;
% end
% xlim([1 100]);
% ylim([0 100]);
% xlabel('# Simulations');
% ylabel('Reflectivity (%)');

% % % Figure 4(d)
% figure; 
% hs=histogram(mValsS,'FaceColor',0.8*[1 1 1],BinWidth=5); 
% hold on; 
% h=histogram(mVals, BinWidth=5);
% xlim([0 100])
% ylim([0 20]);
% xlabel('Reflectivity (%)')
% ylabel('Optimal designs')
% hold on; 
% plot(96.9,1,'b*','MarkerSize',10)

% % % Compare the spacings of a topology-optimization design process at
% % % various iterations
% figure;
% rng(15); % generate same random numbers (to reproduce figure)
% j = find(mVals>97.5);
% ind = j(1);
% figure;
% plot(M{ii}.xHist(2,:));
% hold on;
% plot(M{ii}.xHist(12,:));
% plot(M{ii}.xHist(22,:));
% plot(M{ii}.xHist(32,:));
% plot(M{ii}.xHist(42,:));
