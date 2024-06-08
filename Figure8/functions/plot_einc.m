subplot 231
imagesc(real(einc)); 
title('Re(E{inc})','fontweight','normal');
colorbar('southoutside');
axis equal
xlim([1,Kx])
ylim([1,Ky])

subplot 234
imagesc(imag(einc)); 
title('Im(E{inc})','fontweight','normal');
colorbar('southoutside');
axis equal
xlim([1,Kx])
ylim([1,Ky])

subplot 232
plot(1:Kx,real(einc(sy,:)));hold on;
plot(1:Kx,imag(einc(sy,:)));
plot(1:Kx,real(-1i*besselh(0,1,sqrt(const)*omega*abs(x-x(sx)))/4),'--');
plot(1:Kx,imag(-1i*besselh(0,1,sqrt(const)*omega*abs(x-x(sx)))/4),'--');

plot(PMLx_num*[1,1],[min(real(einc(sy,:))),max(real(einc(sy,:)))],'k:')
plot((Kx-PMLx_num+1)*[1,1],[min(real(einc(sy,:))),max(real(einc(sy,:)))],'k:')

xlabel('x');title('E{inc}','fontweight','normal');
xlim([1,Kx])
ylim(inf*[-1,1])
lgd = legend('Re','Im','location','best');
lgd.ItemTokenSize = [10,18];

subplot 235
plot(1:Ky,real(einc(:,sx)));hold on;
plot(1:Ky,imag(einc(:,sx)));
plot(1:Ky,real(-1i*besselh(0,1,sqrt(const)*omega*abs(y-y(sy)))/4),'--');
plot(1:Ky,imag(-1i*besselh(0,1,sqrt(const)*omega*abs(y-y(sy)))/4),'--');

plot(PMLy_num*[1,1],[min(real(einc(:,sx))),max(real(einc(:,sx)))],'k:')
plot((Ky-PMLy_num+1)*[1,1],[min(real(einc(:,sx))),max(real(einc(:,sx)))],'k:')

xlabel('y');title('E{inc}','fontweight','normal');
xlim([1,Ky])
ylim(inf*[-1,1])
lgd = legend('Re','Im','location','best');
lgd.ItemTokenSize = [10,18];


subplot 233
imagesc(real(es_unstruct));
title('Re(Es^{unstruc.})','fontweight','normal');colorbar('southoutside');

axis equal
xlim([1,Kx])
ylim([1,Ky])

subplot 236
imagesc(imag(es_unstruct));
title('Im(Es^{unstruc.})','fontweight','normal');
colorbar('southoutside');

axis equal
xlim([1,Kx])
ylim([1,Ky])

colormap(jet)
set(findobj(gcf,'type','axes'),'FontName','Calibri',...
    'FontSize',fontSize,'LineWidth',lnWdth,'BoxStyle','full','Box','on');