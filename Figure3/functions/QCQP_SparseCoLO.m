function [Obj, x] = QCQP_SparseCoLO(S0, b0, c0, S, b, c, solvername)
% function [Obj] = QCQP_SparseCoLO(S0, b0, c0, S, b, c, solvername)

global Time

tstart = tic;
disp('Sparse init. start')

m = length(b);
n = length(b0);

% convert QCQP to a SDP
% [F0, F, t] = QCQP_2_SDP(S0, b0, c0, S, b, c);
% 
% n = length(b0);
% f0 = reshape(F0, (n+1)^2, 1);
% f = cell(1, m+1);
% nz = zeros(m+1, 1);
% for i = 1:m+1
%     f{i} = reshape(F{i}, 1, (n+1)^2);
%     F{i} = 0;
%     nz(i) = nnz(f{i});
% end

[f0, f, t] = QCQP_2_SDP(S0, b0, c0, S, b, c);
n = length(b0);
f0 = sparse(reshape(f0, (n+1)^2, 1));
nz = zeros(m+1, 1);
for i = 1:m+1
    f{i} = sparse(reshape(f{i},(n+1)^2, 1));
    nz(i) = nnz(f{i});
end

% --------- experiement ---------
% version 1
%     tic
%     A_old = sparse(m+1, (n+1)^2);
%     for i = 1:m+1
%         A_old(i,:) = f{i}; % this line needs further acceleration
%     end
%     toc

% version 2



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [idx_row, idx_col, idx_val] = deal(zeros(sum(nz), 1));
% idx = 1;
% for i = 1:m+1
%     idx_i = (idx:idx + nz(i)-1);
%     [idx_row_i, idx_col_i, idx_val_i] = find(f{i});
%     idx_row(idx_i) = idx_row_i + i - 1;
%     idx_col(idx_i) = idx_col_i;
%     idx_val(idx_i) = idx_val_i;
%     idx = idx + nz(i);
% end
% A = sparse(idx_row, idx_col, idx_val, m+1, (n+1)^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


idx_col = cellfun(@(x) find(x),f.','UniformOutput',false); % get nonzero indexes
idx_row = zeros(sum(nz), 1);
idx_val = cell(sum(nz), 1);
idx = 1;
for i = 1:m+1
    idx_i = (idx:idx + nz(i)-1);
    idx_row(idx_i) = ones(nz(i),1) * i ;
    idx_val{i} = f{i}(idx_col{i});
    idx = idx + nz(i);
end
A = sparse(idx_row, cell2mat(idx_col), cell2mat(idx_val), m+1, (n+1)^2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     toc

% --------- experiement ---------
b = sparse(m+1,1); 
b(end) = 1;
c = f0;
K.s = n + 1;
J.f = m + 1;

% SparseCoLO + SeDuMi
%[x, y, infoCoLO, cliqueDomain, cliqueRange, LOP] = sparseCoLO(A, b, c, K, J, parCoLO);

parCoLO.SDPsolver = solvername;
parCoLO.domain = 1; parCoLO.range = 2; parCoLO.EQorLMI = 1;

Time.SparseCoLO_init = toc(tstart);
disp('Sparse solve start')

global SG
if SG
    [x, y, infoCoLO, cliqueDomain, cliqueRange, LOP] = sparseCoLO_Mosek_SG(A, b, c, K, J, parCoLO);
    disp('SG on')
else
    [x, y, infoCoLO, cliqueDomain, cliqueRange, LOP] = sparseCoLO_Mosek(A, b, c, K, J, parCoLO);
    disp('SG off')
end




global CLIQUES
CLIQUES.cliqueDomain = cliqueDomain;


[primalObjValue, dualObjValue, primalfeasiblity, dualfeasibility] = ...
    evaluateCoLO(x,y,A,b,c,K,J,cliqueDomain,cliqueRange);
Obj = primalObjValue;

% x = psdCompletion(x, K, cliqueDomain);
x = psdCompletion_SG(x, K, cliqueDomain);

Time.SparseCoLO = toc(tstart) - Time.SparseCoLO_init;
disp('Sparse solve done')

end
