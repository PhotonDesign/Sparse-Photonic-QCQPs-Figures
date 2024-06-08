startup

%% user-defined parameter
w = 1;
chi = (2.3+0.03i)^2-1;
L = 0.3*2*pi;
Nx = 100;
obj = 'R'; % option: R, Pext

%Target phase of reflected wave
targetPhase = pi;

%% init
k = w;

% material
xi = - 1 / chi;

% coordinate
dx = L / (Nx-4);
x = linspace(-2*dx, L+dx, Nx).';
dx = x(2) - x(1);

% Maxwell operator (with absorbing boundary conditions)
M = get_M_ABC(x, w);

% source
rho = sparse(2, 1, 1/dx, Nx, 1);

% calculate incident field 
einc = M \ (- rho);

%% init: objective & constraints
% bound objective 
switch obj
    case 'R'
        D3 = sparse(3,3,1,Nx,Nx);
        S0 = 1/abs(einc(3))^2 * D3;
        b0 = sparse(3,1,1,Nx, 1)*(-2*exp(i*targetPhase))/abs(einc(3));
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

%% cal: bound 
% SparseCoLO + Mosek 
[fobj_bd, x_opt] = complexQCQP_SparseCoLO(S0, b0, c0, Sr, br, cr, 'mosek', 'min');

% extract optimal design
[es_opt, d] = extract_opt(x_opt, Nx);
e_opt = es_opt + einc;
p_opt = - M * es_opt;
chi_opt = p_opt ./ e_opt;

%% reference: thin film
% thin film material
chi_tf = chi * ones(Nx, 1); 
chi_tf([1,2, Nx]) = 0;

e_tf = (M + diag(chi_tf)) \ (- rho);
escat_tf = e_tf - einc;

% reflection
r = escat_tf(3) / einc(3); 
R = abs(r)^2;

% analytical thin film solution
[~, ~, ~, r, ~] = Qana_Q_film(chi + 1, 0, k*L);
R_ana = abs(r)^2;

%% result
fobj_bd
R
R_ana

%% check
NS = length(S);
err = nan(NS, 1);
err_rel = nan(NS, 1);
for i = 1:NS
    err(i) = es_opt' * S{i} * es_opt + b{i}' * es_opt + c{i};
    err_rel(i) = abs(err(i)) ./ max(abs(es_opt' * S{i} * es_opt),  abs(b{i}' * es_opt));
end
err
err_rel

%% test the extracted-optimal design
% structure 1: optimal design
e1 = (M + diag(chi_opt)) \ (- rho);
escat1 = e1 - einc;
r = escat1(3) / einc(3); 
R1 = abs(r)^2;

% structure 2: a more realistic (grey-scale) structure
chi_opt2 = real(chi_opt);
chi_opt2([1,2,Nx]) = 0;
e2 = (M + diag(chi_opt2)) \ (- rho);
escat2 = e2 - einc;
r = escat2(3) / einc(3); 
R2 = abs(r)^2;

% structure 3: binarization
chi_opt3 = zeros(Nx, 1);
chi_opt3(chi_opt2 >= chi / 2) = 1;
e3 = (M + diag(chi_opt3)) \ (- rho);
escat3 = e3 - einc;
r = escat3(3) / einc(3); 
R3 = abs(r)^2;

% plot
figure
set(gcf,'position',[100, 100, 1000, 300])

subplot(1,5,1)
hold on
plot(x, real(chi_tf), 'o-')
plot(x, imag(chi_tf), 'o-')
title(tName({'R'},{R},''))
xlabel('x')
ylabel('\chi')
legend('Re\chi', 'Im\chi','location','south')

subplot(1,5,2)
hold on
plot(x, real(chi_opt), 'o-')
plot(x, imag(chi_opt), 'o-')
title(tName({'R'},{R1},''))

subplot(1,5,3)
hold on
plot(x, real(chi_opt2), 'o-')
plot(x, imag(chi_opt2), 'o-')
title(tName({'R'},{R2},''))

subplot(1,5,4)
hold on
plot(x, real(chi_opt3), 'o-')
plot(x, imag(chi_opt3), 'o-')
title(tName({'R'},{R3},''))

subplot(1,5,5)
hold on
plot(1:length(d), d/d(1), 'o-')
ylabel('eigs(X_{opt}) / \lambda_1')
set(gca, 'yscale','log')

sgtitle(tName({'\omega','R_{tf}', 'R_{bd}'},{w,R_ana, fobj_bd},'N_x = 20 L + 4, '))
