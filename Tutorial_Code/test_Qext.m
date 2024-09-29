% Maximize the extinction coefficient of a multi-layered film
% required: CVX, https://cvxr.com/cvx/download/
%
% Zeyu Kuang (zeyu.kuang@gmail.com) 
% 2024/9/13
startup

%% user-defined parameter
w = 1;
chi = 1; 
L = 10;
Nx = L*50;

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
I = sparse(3:Nx-1, 3:Nx-1, ones(Nx-3,1), Nx, Nx);
S0 = sparse(Nx,Nx);
b0 = -1j*k/abs(einc(3))^2 * (M' * I' * einc) * dx;
c0 = 0;

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
    
%% cal: bound 
Qext_bd = QCQP_solver_complex(S0,b0,c0,S,b,c);

%% result
Qext_bd % theoretical maximum is 4


