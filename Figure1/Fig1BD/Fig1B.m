clear

%% user-defined parameter
d = 0.5;
x1 = [-0.6; -0.5];
x2 = [0.5; 0.2];
x3 = [-0.5; 0.6];
x4 = [0.9; -0.7];
x5 = [-1.2; 1];

fun = @(x) 1 ./ (d + norm(x - x1)) + 1.2 ./ (d + norm(x - x2)) + 1 ./ (d + norm(x - x3)) + 1 ./ (d + norm(x - x4)) + 1 ./ (d + norm(x - x5));

%% cal
x1 = -1:0.01:1;
x2 = -1:0.01:1;

N = length(x1);
y = zeros(N, N);
for i = 1:N
    for j = 1:N
        x = [x1(i); x2(j)];
        y(i, j) = fun(x);
    end
end

%% result
figure
contourf(x1, x2, y, 100,'LineColor','none')

axis square
xticklabels('')
xticks([])
yticklabels('')
yticks([])

set(gcf,'color','w');

%% save
saveas(gcf,'Fig1B.png')
export_fig Fig1B.pdf



