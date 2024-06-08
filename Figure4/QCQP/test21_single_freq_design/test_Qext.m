startup

%% user-defined parameter
w = 1;
chi = 1;
L = 10;
Nx = L*50;
obj = 'R'; % option: R, Qext

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
        b0 = sparse(Nx, 1);
        c0 = 0;
    case 'Qext'
        I = sparse(3:Nx-1, 3:Nx-1, ones(Nx-3,1), Nx, Nx);
        S0 = sparse(Nx,Nx);
        b0 = -1j*k/abs(einc(3))^2 * (M' * I' * einc) * dx;
        c0 = 0;
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
[fobj_bd, x_opt] = complexQCQP_SparseCoLO(S0, b0, c0, Sr, br, cr, 'mosek', 'max');
fobj_bd = real(fobj_bd);

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

% extinction
Qext = real(b0'*escat_tf);

% analytical thin film solution
Qext_ana = Qana_Q_film(chi + 1, 0, k*L);

%% result
fobj_bd
Qext
Qext_ana

%% check
NS = length(S);
err = nan(NS, 1);
for i = 1:NS
    err(i) = es_opt' * S{i} * es_opt + b{i}' * es_opt + c{i};
end
err

%% test the extracted-optimal design
% structure 1: optimal design
e1 = (M + diag(chi_opt)) \ (- rho);
escat1 = e1 - einc;
Qext1 = real(b0'*escat1);

% structure 2: a more realistic (grey-scale) structure
chi_opt2 = real(chi_opt);
chi_opt2([1,2,Nx]) = 0;
e2 = (M + diag(chi_opt2)) \ (- rho);
escat2 = e2 - einc;
Qext2 = real(b0'*escat2);

% structure 3: binarization
chi_opt3 = zeros(Nx, 1);
chi_opt3(chi_opt2 >= chi / 2) = 1;
e3 = (M + diag(chi_opt3)) \ (- rho);
escat3 = e3 - einc;
Qext3 = real(b0'*escat3);

% plot
figure
set(gcf,'position',[100, 100, 1000, 300])

subplot(1,5,1)
hold on
plot(x, real(chi_tf), 'o-')
plot(x, imag(chi_tf), 'o-')
title(tName({'Qext'},{Qext},''))
xlabel('x')
ylabel('\chi')
legend('Re\chi', 'Im\chi','location','south')

subplot(1,5,2)
hold on
plot(x, real(chi_opt), 'o-')
plot(x, imag(chi_opt), 'o-')
title(tName({'Q_{ext}'},{Qext1},''))

subplot(1,5,3)
hold on
plot(x, real(chi_opt2), 'o-')
plot(x, imag(chi_opt2), 'o-')
title(tName({'Q_{ext}'},{Qext2},''))

subplot(1,5,4)
hold on
plot(x, real(chi_opt3), 'o-')
plot(x, imag(chi_opt3), 'o-')
title(tName({'Q_{ext}'},{Qext3},''))

subplot(1,5,5)
hold on
plot(1:length(d), d/d(1), 'o-')
ylabel('eigs(X_{opt}) / \lambda_1')
set(gca, 'yscale','log')

sgtitle(tName({'\omega','Qext_{tf}', 'Qext_{bd}'},{w,Qext_ana, fobj_bd},'N_x = 50 L, '))
