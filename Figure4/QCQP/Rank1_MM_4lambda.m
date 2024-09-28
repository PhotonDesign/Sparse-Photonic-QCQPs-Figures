% pre-requisite: 
% 1. SparseCoLO (http://www.opt.c.titech.ac.jp/kojima/SparseCoLO/SparseCoLO.htm)
% 2. Mosek (https://www.mosek.com/)
%
% After installing SparseCoLO, make sure to add in SparseCoLO112/psdCompletion
% "sDim = K.s(kk);" above line 213.  
%
% Zeyu Kuang
% 2024/6/7

%% Initialize
clc;
clear;
addpath('./1D_VIE_functions')
addpath('./redblue')
addpath('./test21_single_freq_design/function')
addpath('./diff_functions')
%% User-Specified Parameters
n = 2.3+0.03i;
chi = n^2 - 1;
L = 4; % # of wavelengths
Nx = 100; % work: 80, 74. Does not work: 60, 70, 72 
targetPhase = 0.3*pi;

%Modifications for differential QCQP
L_VIE = L; %will need later
L_dif = L*2*pi; %w and k = 1 in differential, not wavelength
phase_dif = pi/2 - targetPhase; %Differential starts out of phase w.r.t VIE

%% First Run
%%% Kept Zeyu's labels as are
%%% user-defined parameter
w = 1;
obj = 'R'; % option: R, Pext

%%% init
k = w;

% material
xi = - 1 / chi;

% coordinate
dx = L_dif / (Nx-4);
x = linspace(-2*dx, L_dif+dx, Nx).';
dx = x(2) - x(1);

% Maxwell operator (with absorbing boundary conditions)
M = get_M_ABC(x, w);

% source
rho = sparse(2, 1, 1/dx, Nx, 1);

% calculate incident field 
einc = M \ (- rho);

%%% init: objective & constraints

% bound objective 
switch obj
    case 'R'
        D3 = sparse(3,3,1,Nx,Nx);
        S0 = 1/abs(einc(3))^2 * D3;
        b0 = sparse(3,1,1,Nx, 1)*(-2*exp(1i*phase_dif))/abs(einc(3));
        c0 = 1;
end

% bound constraints
S = {};
b = {};
c = {};
for i = 1:Nx
    if ismember(i, [1,2,Nx])
        S{end+1} = sparse(Nx, Nx);
        b{end+1} = M(i, :)';
        c{end+1} = sparse(0);
        
        a = M(i,:)';
        S{end+1} = a * a';
        b{end+1} = sparse(Nx,1);
        c{end+1} = sparse(0);
    else
        Di = sparse(i,i,1,Nx,Nx);
        S{end+1} = Di * M - xi' * M' * Di * M;
        b{end+1} = M' * Di' * einc;
        c{end+1} = sparse(0);
    end
end
    
% complex to real constraints
[Sr, br, cr] = c2r_con(S, b, c); 

%%% cal: bound 
% SparseCoLO + Mosek 
[fobj_bd, x_opt, F, F0] = complexQCQP_SparseCoLO(S0, b0, c0, Sr, br, cr, 'mosek', 'min');

% extract optimal design
[es_opt, d,X_opt] = extract_opt(x_opt, Nx);
X_diag = diag(X_opt);
e_opt = es_opt + einc;
p_opt = - M * es_opt;
chi_opt = p_opt ./ e_opt;

%%% test the extracted-optimal design
% structure 1: optimal design
e1 = (M + diag(chi_opt)) \ (- rho);
escat1 = e1 - einc;
r = escat1(3) / einc(3); 
R1 = abs(r)^2;

R0 = (X_diag(3)+X_diag(Nx+3))/abs(einc(3))^2;

%%% Colin's addition
dif_phase = angle(einc(3)) - angle(escat1(3));
if dif_phase > pi
    dif_phase = dif_phase - 2*pi;
elseif dif_phase < -pi
    dif_phase = dif_phase + 2*pi;
end

% Optimization Value
disp(['Optimal reflectivity: ' num2str(R0)])

%Simulated power and phase with first eigenvector
disp(['Simulated reflectivity: ' num2str(R1)])
disp(['Simulated reflection phase: ' num2str(dif_phase/pi) 'pi'])

%Solution Rank
disp(['Ratio of 1st 2 eigenvalues: ' num2str(d(1)/d(2))])

%chiPlot(L_VIE,chi_opt)

[sparse_i,sparse_j,F_graph] = sparsity(F);
dim = length(F{1});

save('Graph_functions/sparsity_test.mat','F_graph')

f0 = reshape(F0,length(F0)^2,1); %Will need later

%% Rank Minimization through MM method
tic
%%%Just do this as one run first to test, then embed in while loop
%%%Stopping condition: eigenvalue ratio of W_opt
eta = 1e-5;
beta = 2; %Factor increase of eta

%Dummies to start
eigRatio = 1;
eigRatio_last = 1e-10;
W_opt_last = 1e-10;

while eigRatio < 1e7

    [A_tr,b_tr,c_tr,K_tr,J_tr] = initialize_trace_pen_2(F0,F,eta);

    x_pen = tr_sparseCoLO_2(A_tr, b_tr, c_tr, K_tr, J_tr);

    [~, d_tr,W_l] = extract_opt(x_pen, Nx);

    %d_tr(1)/d_tr(2)

    del1 = 1e-3;
    del2 = 1e-3;

    d1 = del1 + 1;
    eps = 0.5;
    alpha = 2;
    l = 0;
    while d1 > del1
        k = 0;
        W_k = W_l;
        d2 = del2 + 1;
        while d2 > del2
            W_k_last = W_k;
            g = g_MM_term(W_k_last,eps,eta);
            g = reshape(g,length(g)^2,1);
            c_new = sparse(f0 + g);
            W_k = tr_sparseCoLO_2(A_tr, b_tr, c_new, K_tr, J_tr);
            W_k = reshape(W_k,size(F0));
            d2 = norm(W_k - W_k_last,'fro')/norm(W_k_last,'fro');
            k = k+1;
        end
        W_l_last = W_l;
        W_l = W_k;
        d1 = norm(W_l - W_l_last,'fro')/norm(W_l_last,'fro');
        l = l+1;
        eps = eps/alpha;
    end
    D_l = eigs(W_l);
    eigRatio = D_l(1)/D_l(2);
   %Check if eigenvalue ratio begins to decrease more than 10%
    if eigRatio < eigRatio_last
        W_l = W_opt_last;
        break
    end
    eta = beta*eta
    W_opt_last = W_l;
    eigRatio_last = eigRatio;
end
toc
%% Design from Minimization

x_opt = reshape(W_l(1:(2*Nx+1),1:(2*Nx+1)),(2*Nx+1)^2,1);
[es_opt, d,X_opt] = extract_opt(x_opt, Nx);
X_diag = diag(X_opt);
e_opt = es_opt + einc;
p_opt = - M * es_opt;
chi_opt = p_opt ./ e_opt;

%%% test the extracted-optimal design
% structure 1: optimal design
e1 = (M + diag(chi_opt)) \ (- rho);
escat1 = e1 - einc;
r = escat1(3) / einc(3); 
R1 = abs(r)^2;

R0 = (X_diag(3)+X_diag(Nx+3))/abs(einc(3))^2;

%%% Colin's addition
dif_phase = angle(einc(3)) - angle(escat1(3));
if dif_phase > pi
    dif_phase = dif_phase - 2*pi;
elseif dif_phase < -pi
    dif_phase = dif_phase + 2*pi;
end

% Optimization Value
disp(['Optimal reflectivity: ' num2str(R0)])

%Simulated power and phase with first eigenvector
disp(['Simulated reflectivity: ' num2str(R1)])
disp(['Simulated reflection phase: ' num2str(dif_phase/pi) 'pi'])

%Solution Rank
disp(['Ratio of 1st 2 eigenvalues: ' num2str(d(1)/d(2))])

%chiPlot_pos(L_VIE,chi_opt) 

%% Perfect Binarization
chi_threshold = 0.4;
valType = 'abs';

chi_dif_b = binaryDesign(chi,chi_opt,chi_threshold,valType,Nx);
[R_b,~,~,dif_phase_b] = getPrefFromChi(chi_dif_b,L,Nx,1);

disp(['Binary R: ' num2str(R_b)])
disp(['Binary phase: ' num2str(dif_phase_b/pi) 'pi'])

%chiPlot_b(L,chi_dif_b)

%% Result (ZK)
lam = 2*pi;

figure
set(gcf, 'position', [100, 100, 1000, 400])

subplot(1,2,1)
hold on
plot(x / lam, real(chi_opt))
plot(x / lam, imag(chi_opt))
xlabel('x / \lambda')
ylabel('\chi')
legend('real', 'imag')
axis tight


subplot(1,2,2)
hold on
plot(x / lam, real(chi_dif_b))
plot(x / lam, imag(chi_dif_b))
xlabel('x / \lambda')
ylabel('\chi')
legend('real', 'imag')
title('binarized')
axis tight

%% save
%save('Rank1_MM')