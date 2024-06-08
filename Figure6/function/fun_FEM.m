
function [xx, e_FEM] = fun_FEM(N, Nx)
    % N: resolution of FEM
    % Nx: resolution of plot

    rng(13)
    
    %% user-defined parameter
    k = 1;
%     N = 10000;
%     Nx = 15731;
    L = 8*pi+0.4;
    s = 2; % source location
    chi_mat = 11;
    Nm = 10; % number of material patches
    
    %% init
    xm = sort(rand(2*Nm,1)) * L;
    xm(1) = 0.2;
    xm(end) = xm(1)+8*pi;
    xm(end-5) = 19.2;
    xm = xm - 0.2;
    Xm = [xm(1:2:end), xm(2:2:end)]; % coordinates for the material patches
    
    x = linspace(-0.2, L-0.2, N).';
    x = x-min(abs(x));
    h = x(2) - x(1);
    
    gamma1 = -1j*k; % boundary terms
    gamma2 = 1j*k;
    
    chi = zeros(N, 1);
    for i = 1:N
        if any(x(i) >= Xm(:,1) & x(i) < Xm(:,2)) % the point is inside a material patch
            chi(i) = chi_mat * any(x(i) >= Xm(:,1) & x(i) < Xm(:,2));
        end
    end 
    
    %% cal
    K1 = 1 / h * sparse([1:N-1, 2:N, 1:N], [2:N, 1:N-1, 1:N], [-1*ones(1, N-1), -1*ones(1, N-1), 1, 2*ones(1,N-2), 1], N, N);
    K2 = k^2 * h / 6 * sparse([1:N-1, 2:N, 1:N], [2:N, 1:N-1, 1:N], [ones(1, N-1), ones(1, N-1), 2, 4*ones(1,N-2), 2], N, N);
    K4 = sparse([1, N], [1, N], [gamma1, -gamma2], N, N);
    
    [left_edge_pos, right_edge_pos] = transform_chi(chi);
    K3 = generate_K3(left_edge_pos, right_edge_pos, chi_mat, h, N);
    
    K = K1 - K2 - K3 + K4;
    b = 1j*k*sparse(s, 1, 1, N, 1).* exp(-1j*((0-x(s))*k+pi));
    
    
    u = K\b;
    
    %% post processing
    % FEM solution
    tot_thick = xm(end)-xm(1);
    xx = linspace(0,tot_thick,Nx);
    Phi = zeros(Nx, N);
    phi1 = @(x) (1 - x/h) .* (0 < x & x <= h);
    phi2 = @(x) (1 + x/h) .* (-h <= x & x <= 0);
    for i = 1:N
        if i == 1
            Phi(:, i) = phi1(xx - x(1) +eps);
        elseif i == N
            Phi(:, i) = phi2(xx - x(N));
        else
            Phi(:, i) = phi1(xx - x(i)) + phi2(xx - x(i));
        end
    end
    
    e_FEM = Phi * u;
    
    %% function
    
    function [left_edge_pos, right_edge_pos] = transform_chi(chi)
        Ntemp = length(chi);
        left_edge_pos = [];
        right_edge_pos = [];
        for ii = 1:Ntemp-1
            if chi(ii)==0 && chi(ii+1)~=0
                left_edge_pos = [left_edge_pos; ii+1];
            elseif chi(ii)~=0 && chi(ii+1)==0
                right_edge_pos = [right_edge_pos; ii+1];
    
            end
        end
    end
    
    function K3 = generate_K3(left_edge_pos, right_edge_pos, chi_mat, h, N)
        temp = chi_mat*h*[1/3, 1/6; 1/6, 1/3];
        K3 = zeros(N,N);
        N_edge = length(left_edge_pos);
        for ii = 1:N_edge
            for j = left_edge_pos(ii):right_edge_pos(ii)-1
                K3(j:j+1, j:j+1) = K3(j:j+1, j:j+1) + temp;
            end
        end
    end

end











