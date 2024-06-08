
function[] = createMovie(xHist,mHist,fc,Niter)

xmax = 1.02*max(sum(xHist,2));
xv = linspace(0,xmax,300);

for i = 1:Niter
epsr = plotDesign(xHist(i,:),xv,3,0);
figure(1);
sbp1 = subplot(2,1,1);
plot(fc(1:i),mHist(1:i),'k');
xlim([1 max(max(fc),2)]);
ylim([0 1]);
yticks(0:0.1:1);
yticklabels({'0','','','','','0.5','','','','','1'})
xlabel('# Simulations');
ylabel('Reflectivity');
sbp2 = subplot(2,1,2);
imagesc(epsr)
colormap gray
axis off
sbp2.Position = [0.13,0.136,0.775,0.1];
sbp1.Position = [0.13,0.4,0.775,0.5];
F(i) = getframe(gcf) ;
drawnow
end

writerObj = VideoWriter('myVideo.mp4', 'MPEG-4');
writerObj.FrameRate = 10;
open(writerObj);
for i=1:length(F)   
    writeVideo(writerObj, F(i));
end
close(writerObj);
end