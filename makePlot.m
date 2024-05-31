close all; 

%% load data
fig = openfig('Sweep length2.fig');

dataObjs = findobj(fig,'-property','YData');
y1 = dataObjs(3).YData;
y2 = dataObjs(4).YData;
y3 = dataObjs(5).YData;
y4 = dataObjs(6).YData;
x = dataObjs(3).XData;


figure
plot(x,y1); hold on
plot(x,y2)
plot(x,y3)
plot(x,y4)


%% new plot
set(0,'DefaultLineLineWidth', 1);
set(0,'DefaultAxesTitleFontWeight','normal');

% colors
c1 = [228,26,28]./255;
c2 = [55,126,184]./255;
c3 = [77,175,74]./255;
c4 = [152,78,163]./255;

c5 = [255,127,0]./255;
c6 = [255,255,51]./255;
c7 = [166,86,40]./255;
c8 = [247,129,191]./255;

figure(10);set(gcf,'Position',[800,600,300,250]);clf;hold on

n = 1.1153e+06;

ll = 1:(23);
% xx = x(ll)./x(ll(end));
xx = x(ll);

plot(xx,y1(ll)./n,'o-','Color',c1)
plot(xx,y2(ll)./n,'d-','Color',c2)
plot(xx,y3(ll)./n,'s-','Color',c3)
plot(xx,y4(ll)./n,'x-','Color',c4)

lgd=legend('FEM bound','FDFD bound','FEM sim.','FDFD sim.','location','best');
lgd.ItemTokenSize = [12,18];

xlabel('{\it{L}}_{des}  (\mum)');
ylabel('{\it f} ({\bfE})  (a.u.)')

xlim([-0.01,.93])
ylim([-.02,1.075])

set(findobj(gcf,'type','axes'),'FontName','Calibri',...
    'FontSize',14,'LineWidth',1,'BoxStyle','full','Box','on');


