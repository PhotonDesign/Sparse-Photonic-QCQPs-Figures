function[layerMaterials] = createLayerMaterials(eps_array, layerThicknesses)
layerMaterials = cell(1,length(eps_array)+2);
layerMaterials{1} = @Vacuum;
layerMaterials{end} = @Vacuum;
for i=1:length(layerThicknesses)
    layerMaterials{i+1} = @(w)eps_array(i);
end
end