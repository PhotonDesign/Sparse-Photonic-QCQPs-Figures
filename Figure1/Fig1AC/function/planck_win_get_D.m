function [D0,D1] = planck_win_get_D(kin,hin)


h = 1/8;
k = kin*hin/h;
a = 10;
N = a/h;


err = 0.49;
wl = @(n)(1./(1+exp(err*2*N./n - err*2*N./(err*2*N-n)))).*(n>=1).*(n<err*2*N) + (n<=N).*(n>=err*2*N);
ws = @(n) wl(n) + wl(2*N-n).*(n>(2*N/2));
w = @(n) ws(n+N);

%wind = @(x,y) w(x/h).*w(y/h);
wind = @(x,y) w(sqrt((x/h).^2 + (y/h).^2));
%wind = @(x,y) exp(-(x.^2 + y.^2)*40);

fD0 = @(x,y) besselh(0, k*sqrt(x.^2 + y.^2));


wfD0 = @(x,y) fD0(x,y).*wind(x,y);
wfD1 = @(x,y) (x.^2).*besselh(0, k*sqrt(x.^2 + y.^2)).*wind(x,y);
wfD2 = @(x,y) (x.^4).*besselh(0, k*sqrt(x.^2 + y.^2)).*wind(x,y);
x = -a:h:a;y = -a:h:a;

% figure;
% plot(x,wind(x,0));
% xx = linspace(-h,h,100);
% figure;
% plot(xx,wind(xx,0));

%%

temp = ones(length(x),1); 
wx = temp*h;
temp = ones(length(y),1); 
wy = temp*h;
W = (wy).*(wx.');

digits(32); 
[X2d,Y2d] = meshgrid(x,y); hval = wfD0(X2d,Y2d); hval(N+1,N+1) = 0;
inttrap = sum(W.*hval, 'all');

intex = integral2(wfD0,0,a,0,a,'RelTol', 1e-16);
D0 = (4*intex - inttrap)/(h^2);
%disp('h = '); disp(h);
%disp('D0 = '); disp(D0);

[X2d,Y2d] = meshgrid(x,y); hval = wfD1(X2d,Y2d); hval(N+1,N+1) = 0;
inttrap = sum(W.*hval, 'all');

intex = integral2(wfD1,0,a,0,a,'RelTol', 1e-16);
%intex = integral2(wfD1,0,a,0,a);
D1 = (4*intex - inttrap)/(h^4);
%disp('h = '); disp(h);
%disp('D1 = '); disp(D1);

% [X2d,Y2d] = meshgrid(x,y); hval = wfD2(X2d,Y2d); hval(N+1,N+1) = 0;
% inttrap = sum(W.*hval, 'all');
% 
% intex = integral2(wfD2,0,a,0,a,'RelTol', 1e-16);
% %intex = integral2(wfD1,0,a,0,a);
% D2 = (4*intex - inttrap)/(h^6);
% disp('h = '); disp(h);
% disp('D2 = '); disp(D2);


digits(16); 
end

