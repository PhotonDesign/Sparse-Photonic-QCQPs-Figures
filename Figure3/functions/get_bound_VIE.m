function[f_opt,X] = get_bound_VIE(D,AA,BB,quiet)

Nx = size(AA,1)-1;

if quiet
    cvx_begin quiet
else
    cvx_begin
end

% cvx_solver mosek

variable X(Nx+1,Nx+1) hermitian

AA_trans = AA.';
maximize( real(   AA_trans(:).'*X(:)   ) ); % equivalent to: maximize(trace(AA*P));

subject to

for ii = 1:sum(~cellfun('isempty',D))
    tmp = (BB(D{ii})).';
    tmp(:).'*X(:) == 0; % equivalent to: trace(BB(D{ii})*P) == 0   
end

X == hermitian_semidefinite(Nx+1);
X(Nx+1,Nx+1) == 1;

cvx_end

f_opt = cvx_optval;

end