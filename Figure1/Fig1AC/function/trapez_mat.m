function B = trapez_mat(ax, ay, n, m)
    hx = ax / (n-1);
    hy = ay / (m-1);

    x = 0:hx:ax;
    y = 0:hy:ay;
    [X, Y] = meshgrid(x,y);
    corner_index = (X == 0 & Y == 0) | (X == 0 & Y == ay) ...
                   | (X == ax & Y == 0) | (X == ax & Y == ay);
    edge_index = xor(X == 0 | X == ax | Y == 0 | Y == ay, corner_index);
    B = ones(size(X));
    B(corner_index) = 1/4;
    B(edge_index) = 1/2;
end