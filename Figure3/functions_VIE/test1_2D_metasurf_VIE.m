%% calculate bounds using VIE

AA = [A_VIE,1/2*b_VIE;1/2*b_VIE',c_VIE]; % homogenized matrix with target function

BB = @(D) [Re(D*omega'*(G0_VIE+(xi*eye(Kx_VIE*Ky_VIE)))),1/2*omega'*D*einc_vect_VIE ;...
    1/2*omega*einc_vect_VIE'*D' , 0]; % constraints

D = get_all_D(Kx_VIE*Ky_VIE); % using all D matrices


VIE_solve_tic = tic;
%get bounds
[fobj_bd_VIE,~] = get_bound_VIE(D,AA,BB,1); % last argument sets cvx to quiet






%%
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if VIE_plot_ON
    
    %% check that GF is correct (incident field)
    tmp_idx = reshape((1:Kx*Ky).',[Kx,Ky]).';
    tmp_idx = tmp_idx(sy,sx);
    tmpG0 = reshape(G0_2D(tmp_idx,:),[Kx,Ky]).';
    
    %
    figure(501);set(gcf,'Position',[400,100,300,600]);clf;fontSize = 14;lnWdth = 1;
    
    clim_real_max = max(max(real(einc(:))),max(max(real(tmpG0))));
    clim_real_min = min(min(real(einc(:))),min(min(real(tmpG0))));
    
    subplot 221;
    imagesc(real(einc));title('Re(E{inc}) Diff.','fontweight','normal');
    colorbar('southoutside');axis square;caxis([clim_real_min,clim_real_max])
    
    subplot 222;
    imagesc(real(tmpG0));title('Re(E{inc}) VIE','fontweight','normal');
    colorbar('southoutside');axis square;caxis([clim_real_min,clim_real_max])
    
    subplot 223;
    imagesc(Nx+1:Mx,Ny+1:My,real(einc(Ny+1:My,Nx+1:Mx)));title('Re(E{inc}) Diff.','fontweight','normal');
    colorbar('southoutside');axis square;caxis([clim_real_min,clim_real_max])
    
    subplot 224;
    imagesc(Nx+1:Mx,Ny+1:My,real(tmpG0(Ny+1:My,Nx+1:Mx)));title('Re(E{inc}) VIE','fontweight','normal');
    colorbar('southoutside');axis square;caxis([clim_real_min,clim_real_max])
    
    colormap(jet)
    set(findobj(gcf,'type','axes'),'FontName','Calibri',...
        'FontSize',fontSize,'LineWidth',lnWdth,'BoxStyle','full','Box','on');
    
    %%
    figure(509);set(gcf,'Position',[800,100,300,600]);clf;fontSize = 14;lnWdth = 1;
    
    clim_imag_max = max(max(imag(einc(:))),max(max(imag(tmpG0))));
    clim_imag_min = min(min(imag(einc(:))),min(min(imag(tmpG0))));
    
    subplot 221;
    imagesc(imag(einc));title('Im(E{inc}) Diff.','fontweight','normal');
    colorbar('southoutside');axis square;caxis([clim_imag_min,clim_imag_max])
    
    subplot 222;
    imagesc(imag(tmpG0));title('Im(E{inc}) VIE','fontweight','normal');
    colorbar('southoutside');axis square;caxis([clim_imag_min,clim_imag_max])
    
    subplot 223;
    imagesc(Nx+1:Mx,Ny+1:My,imag(einc(Ny+1:My,Nx+1:Mx)));title('Im(E{inc}) Diff.','fontweight','normal');
    colorbar('southoutside');axis square;caxis([clim_imag_min,clim_imag_max])
    
    subplot 224;
    imagesc(Nx+1:Mx,Ny+1:My,imag(tmpG0(Ny+1:My,Nx+1:Mx)));title('Im(E{inc}) VIE','fontweight','normal');
    colorbar('southoutside');axis square;caxis([clim_imag_min,clim_imag_max])
    
    colormap(jet)
    set(findobj(gcf,'type','axes'),'FontName','Calibri',...
        'FontSize',fontSize,'LineWidth',lnWdth,'BoxStyle','full','Box','on');
    
    %% incident field
    figure(502);set(gcf,'Position',[800,300,300,400]);clf;fontSize = 14;lnWdth = 1;
    
    subplot 221;
    plot(real(tmpG0(sy,:)));hold on
    plot(imag(tmpG0(sy,:)))
    plot(1:Kx,real(einc(sy,:)),'--');
    plot(1:Kx,imag(einc(sy,:)),'--');
    lgd = legend('Re','Im','location','best');lgd.ItemTokenSize = [10,18];
    xlabel('x');xlim(inf*[-1,1]);title('E{inc}','fontweight','normal')
    
    subplot 222;
    plot(real(tmpG0(:,sx)));hold on
    plot(imag(tmpG0(:,sx)));
    plot(1:Ky,real(einc(:,sx)),'--');
    plot(1:Ky,imag(einc(:,sx)),'--');
    lgd = legend('Re','Im','location','best');lgd.ItemTokenSize = [10,18];
    xlabel('y');xlim(inf*[-1,1]);title('E{inc}','fontweight','normal')
    
    subplot 223;
    plot(Nx+1:Mx,real(tmpG0(sy,Nx+1:Mx)));hold on
    plot(Nx+1:Mx,imag(tmpG0(sy,Nx+1:Mx)))
    plot(Nx+1:Mx,real(einc(sy,Nx+1:Mx)),'--');
    plot(Nx+1:Mx,imag(einc(sy,Nx+1:Mx)),'--');
    lgd = legend('Re','Im','location','best');lgd.ItemTokenSize = [10,18];
    xlabel('x');xlim(inf*[-1,1]);title('E{inc}','fontweight','normal')
    
    subplot 224;
    plot(Ny+1:My,real(tmpG0(Ny+1:My,sx)));hold on
    plot(Ny+1:My,imag(tmpG0(Ny+1:My,sx)));
    plot(Ny+1:My,real(einc(Ny+1:My,sx)),'--');
    plot(Ny+1:My,imag(einc(Ny+1:My,sx)),'--');
    lgd = legend('Re','Im','location','best');lgd.ItemTokenSize = [10,18];
    xlabel('y');xlim(inf*[-1,1]);title('E{inc}','fontweight','normal')
    
    set(findobj(gcf,'type','axes'),'FontName','Calibri',...
        'FontSize',fontSize,'LineWidth',lnWdth,'BoxStyle','full','Box','on');
    
    
    
    %% unstructured dielectric with VIE
    figure(503);set(gcf,'Position',[1200,100,300,600]);clf;fontSize = 14;lnWdth = 1;
    
    clim_real_max = max(max(real(es_unstruct_VIE(:))),max(max(real(es_unstruct(Ny+1:My,Nx+1:Mx)))));
    clim_real_min = min(min(real(es_unstruct_VIE(:))),min(min(real(es_unstruct(Ny+1:My,Nx+1:Mx)))));
    clim_imag_max = max(max(imag(es_unstruct_VIE(:))),max(max(imag(es_unstruct(Ny+1:My,Nx+1:Mx)))));
    clim_imag_min = min(min(imag(es_unstruct_VIE(:))),min(min(imag(es_unstruct(Ny+1:My,Nx+1:Mx)))));
    
    subplot 221;
    imagesc(Nx+1:Mx,Ny+1:My,real(es_unstruct(Ny+1:My,Nx+1:Mx)));colorbar('southoutside');
    title('Re(E{s}^{unstruc.}) Diff.','fontweight','normal');axis square
    caxis([clim_real_min,clim_real_max])
    
    subplot 222;
    imagesc(Nx+1:Mx,Ny+1:My,real(es_unstruct_VIE));colorbar('southoutside');
    title('Re(E{s}^{unstruc.}) VIE','fontweight','normal');axis square
    caxis([clim_real_min,clim_real_max])
    
    subplot 223;
    imagesc(Nx+1:Mx,Ny+1:My,imag(es_unstruct(Ny+1:My,Nx+1:Mx)));colorbar('southoutside');
    title('Im(E{s}^{unstruc.}) Diff.','fontweight','normal');axis square
    caxis([clim_imag_min,clim_imag_max])
    
    subplot 224;
    imagesc(Nx+1:Mx,Ny+1:My,imag(es_unstruct_VIE));colorbar('southoutside');
    title('Im(E{s}^{unstruc.}) VIE','fontweight','normal');axis square
    caxis([clim_imag_min,clim_imag_max])
    
    colormap(jet)
    set(findobj(gcf,'type','axes'),'FontName','Calibri',...
        'FontSize',fontSize,'LineWidth',lnWdth,'BoxStyle','full','Box','on');
    
end