% This script: simulate EM scattering with a method based on volume-integral equation.
clear
startup
rng(5)

%% user-defined parameters
% square domain
ax = 0.5*pi; % side length in x
ay = 0.5*pi; % side length in y
nx = 9; % number of points in the x direction
ny = 9; % number of points in the y direction

% material chi
NG = 2; % number of Gaussian bumps
G_max = 25; % Gaussian max
G_width = 0.2; % Gaussian width

% quadrature order
quad_order = 6; % option 0, 4, 6

%% initialization
k = 1;

% coordinates
hx = ax / (nx-1);
hy = ay / (ny-1);
x = 0:hx:ax;
y = 0:hy:ay;
[X, Y] = meshgrid(x,y);
xy = [X(:), Y(:)];
Ngrid = length(xy);
dS = hx * hy;

% susceptibility matrix
G_xy = [rand(NG,1) * ax, rand(NG,1) * ay]; % x-y coodinates of the center of the Gaussian bumps
chi_func = @(x, y) 0;
for i = 1:NG
    chi_func = @(x,y) chi_func(x, y) + exp(-((x - G_xy(i, 1)).^2 + (y - G_xy(i, 2)).^2) / G_width) * G_max; % susceptibility distribution
end

chi_vec = chi_func(xy(:,1), xy(:,2));

% binarization
th = max(chi_vec) / 4;
chi_vec(chi_vec >= th) = max(chi_vec);
chi_vec(chi_vec < th) = 0;

Chi = diag(chi_vec);

% incident field
einc = exp(1j*k*xy(:,1)); % unit plane wave in +x direction

% Green's function (with quadrature correction)
g = GreensFunction_Kernel(k, ax, ay, hx, hy, quad_order);
G = t2BTTB(g); % Green's function matrix

% trapezoidal matrix
bmat = trapez_mat(ax, ay, nx, ny);
b = reshape(bmat, nx*ny, 1);
B = diag(b);

%% calculation: homogeneous material
% direct inverse
I = eye(Ngrid);
p = - (Chi * G * B * dS - I) \ (Chi * einc);

% Extinction coefficient
pftor = k * dS / ay; % pre-fractor 
Qext = pftor * imag(einc' * B * p); 

% fields
escat = G * B * dS * p; % scattered field
e = einc + escat; % total field


%% post-processing: interpolation
% interpolating fields in the volume
xv = 0:hx:ax;
yv = 0:hy:ay;
[Xv, Yv] = meshgrid(xv, yv);
E = real(reshape(e, ny, nx));
Ev = interp2(X,Y,E,Xv,Yv);

%% result
ct = G_max / 2;
cmin = min(real(Ev(:)));
cmax = max(real(Ev(:)));
dx = x(2) - x(1);
xx = linspace(x(1) - dx/2, x(end) + dx/2, nx + 1);

figure 
ax1 = subplot(1,1,1);
hold on
imagesc(x, y, reshape(chi_vec, ny, nx))
colormap(ax1, flipud(gray))
caxis([0,40])
axis off
axis tight
axis equal

hm = mesh(xx, xx, zeros(nx+1, nx+1));
hm.FaceColor = 'none';
hm.EdgeColor = 'k';


saveas(gcf,'Fig1A.png')
export_fig Fig1A.pdf


figure
ax2 = subplot(1,1,1);
hold on
imagesc(xv, yv, real(Ev))
caxis([cmin,cmax])
colormap(ax2, meep)
axis off
axis tight
axis equal

hm = mesh(xx, xx, zeros(nx+1, nx+1));
hm.FaceColor = 'none';
hm.EdgeColor = 'k';

set(gcf, 'color', 'w')

saveas(gcf,'Fig1C.png')
export_fig Fig1C.pdf




