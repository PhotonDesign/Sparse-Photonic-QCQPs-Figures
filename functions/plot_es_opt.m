subplot 211
semilogy(svd(X)/max(svd(X)),'o-','linewidth',lnWdth);axis square
xlim([-.1*length(X),1.1*length(X)]);ylim([-inf,5])
yticks(logspace(-10,0,3))
title('log[ eig(X) ]','fontweight','normal')

subplot 212
imagesc(real(es_opt)); title('Re(Es^{opt}) Diff.','fontweight','normal');axis square

set(findobj(gcf,'type','axes'),'FontName','Calibri',...
    'FontSize',fontSize,'LineWidth',lnWdth,'BoxStyle','full','Box','on');