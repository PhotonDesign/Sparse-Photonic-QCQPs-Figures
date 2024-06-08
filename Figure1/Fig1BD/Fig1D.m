A = [1 0.4;
     0 1];

A = A / norm(A);

x1 = -1:0.01:1;
x2 = -1:0.01:1;

N = length(x1);
y = zeros(N, N);
for i = 1:N
    for j = 1:N
        x = [x1(i); x2(j)];
        y(i, j) = x'*A*x;
    end
end

figure
contourf(x1, x2, y, 20,'LineColor','none')


axis square
xticklabels('')
xticks([])
yticklabels('')
yticks([])

set(gcf,'color','w');
saveas(gcf,'Fig1D.png')
export_fig Fig1D.pdf




