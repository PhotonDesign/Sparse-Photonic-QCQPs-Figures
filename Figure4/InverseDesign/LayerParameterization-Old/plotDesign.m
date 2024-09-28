% lt = layer thicknesses
% xv, Ny
function[epsr] = plotDesign(lt,xv,Ny,plotDesign)
Nt = length(lt);
epsr = ones(size(xv));
x0 = xv(1);
for i=1:Nt
    ind = logical((xv>=x0).*(xv<x0+lt(i)));
    if (mod(i,2)==1)
        epsr(ind) = 0;
    end
    x0 = x0 + lt(i);
end
epsr = repmat(epsr,Ny,1);
if (plotDesign)
    figure; imagesc(epsr)
    colormap gray
    axis off
    set(gcf,'Position',[200 800   670    55])
end
end