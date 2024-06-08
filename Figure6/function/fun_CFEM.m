function [xx, e_FEM] = fun_CFEM(N, Nx, width_factor)
% N: resolution of FEM
% Nx: resolution of plot
% width_factor: width factor of material

rng(13)

%% user-defined parameter
k = 1;
% N = 1000;
% Nx = 15731;
L = 8*pi+0.4;
s = 2; % source location
chi_mat = 11;
% width_factor = [1];
Nm = 10; % number of material patches

%% init
xm = sort(rand(2*Nm,1)) * L;
xm(1) = 0.2;
xm(end) = xm(1)+8*pi;
xm(end-5) = 19.2;
Xm = [xm(1:2:end), xm(2:2:end)]; % coordinates for the material patches

%% cal
Nw = length(width_factor);
[e_FEM, xx] = CFEM(k, N, Nx, L, s, chi_mat, width_factor, Xm);
xx = xx-xx(1);

%% function
function [e_FEM, xx] = CFEM(k, N, Nx, L, s, chi_mat, width_factor, Xm)
    %% init
    x = linspace(0, L, N).';
    h = x(2) - x(1);
    w = width_factor * h / 2; % material width
    
    gamma1 = -1j*k; % boundary terms
    gamma2 = 1j*k;
    
    K1 = 1 / h * sparse([1:N-1, 2:N, 1:N], [2:N, 1:N-1, 1:N], [-1*ones(1, N-1), -1*ones(1, N-1), 1, 2*ones(1,N-2), 1], N, N);
    K2 = k^2 * h / 6 * sparse([1:N-1, 2:N, 1:N], [2:N, 1:N-1, 1:N], [ones(1, N-1), ones(1, N-1), 2, 4*ones(1,N-2), 2], N, N);
    K4 = sparse([1, N], [1, N], [gamma1, -gamma2], N, N);
    
    chi = zeros(N-2, 1);
    for i = 2:N-1
        if any(x(i) > Xm(:,1) & x(i) < Xm(:,2)) % the point is inside a material patch
            chi(i) = chi_mat * any(x(i) > Xm(:,1) & x(i) < Xm(:,2));
        end
    end 
    chi = [0; 0; chi; 0; 0];
    K3 = zeros(N, N); 
    for i = 1:N
        ii = i + 1; 
        K3(i, i) = (chi(ii-1) + chi(ii+1)) * w^2 / 6 / h + chi(ii) * h * (1/3 * w^2 / h^2 - w / h + 1);
        if i <= N-1
            K3(i, i+1) = (chi(ii) + chi(ii+1)) * w / 2 / h * (h/2 - w/3);
            K3(i+1, i) = K3(i, i+1);
        end
    end
    
    K = K1 - K2 - K3 + K4;
    b = 1j*k*sparse(s, 1, 1, N, 1).* exp(-1j*((0.2-x(s))*k+pi));
    
    %% cal
    u = K \ b;
    
    %% post processing
    % FEM solution
%     xx = linspace(0, L, Nx).';
    xx = linspace(0.2,L-0.2,Nx).'-h/2;
    Phi = zeros(Nx, N);
    phi1 = @(x) (1 - x/h) .* (0 < x & x <= h);
    phi2 = @(x) (1 + x/h) .* (-h <= x & x <= 0);
    for i = 1:N
        if i == 1
            Phi(:, i) = phi1(xx+eps);
        elseif i == N
            Phi(:, i) = phi2(xx - x(N));
        else
            Phi(:, i) = phi1(xx - x(i)) + phi2(xx - x(i));
        end
    end
    
    e_FEM = Phi * u;
end 



end