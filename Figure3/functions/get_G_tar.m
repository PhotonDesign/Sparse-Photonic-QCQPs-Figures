function[G_tar] = get_G_tar(omega,x,y,x_tar,y_tar)

[X1,Y1] = meshgrid(x,y);
X1_vect = reshape(X1.',[],1);
Y1_vect = reshape(Y1.',[],1);

r2 = (X1_vect - x_tar).^2 + (Y1_vect - y_tar).^2;

G_tar = (-1i/4).* besselh(0,1,omega*sqrt(r2)).';

% still need to multiply by -k^2*h^2 for our convention
