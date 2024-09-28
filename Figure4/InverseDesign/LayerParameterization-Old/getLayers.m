
function[m] = getLayers(x,x0,w,layerThicknesses,layerMaterials)
m = ones(1,length(x));
N = length(layerThicknesses);
xStart = x0;
for i=1:N
    xEnd = xStart + layerThicknesses(i);
    cond = logical((x>=xStart).*(x<xEnd));
    m(cond) = layerMaterials{i}(w);
    xStart = xEnd;
end
m = [m; m];
end