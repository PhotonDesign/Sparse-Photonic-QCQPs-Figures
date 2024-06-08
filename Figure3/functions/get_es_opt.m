function [es_opt,es_opt_vect] = get_es_opt(X)

[es_opt_tmp,es_v_opt] = eigs(X,1); % estimate optimal p using the eigenvector of largest eigenvalue
es_opt_tmp = sqrt(es_v_opt)*es_opt_tmp;
es_opt_vect = es_opt_tmp(1:end-1)/es_opt_tmp(end); % "undoing" the phase of |s|=1

Nx = sqrt(size(X,1)-1);
My = Nx;

es_opt = reshape(es_opt_vect,[Nx,My]).';