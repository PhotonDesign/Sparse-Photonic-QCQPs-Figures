
function[] = plotLayers(m)
    imagesc(m);
    mm = max(abs(m(:)));
    caxis([-mm,mm]+1);
    colormap(meep);
end