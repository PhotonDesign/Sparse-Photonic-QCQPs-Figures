function Wm = get_COMSOL_tri_trapz_mat(xminfo,COMSOL_idx,units)

if nargin==3
    if strcmp(units,'um')
        scale = 1e6;
    end
else
scale = 1;    
end


% 2D trapezoidal integration matrix using triangles from COMSOL mesh

tris = xminfo.elements.tri.dofs;
coords = xminfo.dofs.coords;
xx = coords(1,:).'*scale;
yy = coords(2,:).'*scale;

NN = size(tris,2);

w = zeros(1,length(xx));

big_tri_corner_idx = [1,2,5;2,3,4;2,4,5;4,5,6];

for pp = 1:NN
    
    for rr = 1:4

        corners_idx = tris(big_tri_corner_idx(rr,:),pp)+1;
        
        if sum(ismember(corners_idx,COMSOL_idx))==3
            area = abs(1/2*det([xx(corners_idx),yy(corners_idx),[1;1;1]]));
            w(corners_idx) = w(corners_idx) + 1/3*area;
        end       
    end
end

Wm = diag(w);

end
