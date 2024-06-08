disp(['L_tot = ',num2str(L),', ',...
    'h = ',num2str(h),', ',...
    'Nx = ',num2str(K),', ',...
    'lambda = ',num2str(lambda),', ',...
    'PML thickness = ',num2str(BC{1}{2}(1)),', '...
    'PML beta = ',num2str(BC{1}{2}(2))])

disp(['L_mat = ',num2str(L_mat),', ',...
    'chi = ',num2str(chi),', ',...
    'obj = ',num2str(obj_no)])

disp(' ')

disp(['cvx + Mosek time: ',num2str(time_cvx),' s']);
disp(['SparseCoLo + Sedumi time: ',num2str(time_SparseCoLoSedumi),' s']);
disp(['SparseCoLo + Mosek time: ',num2str(time_SparseCoLoMosek),' s']);
disp(' ')

disp(['cvx + Mosek bound: ',num2str(round(fobj_bd_cvx,9))])
disp(['SparseCoLo + Sedumi bound: ',num2str(round(fobj_bd_SparseCoLoSedumi,9))])
disp(['SparseCoLo + Mosek bound: ',num2str(round(fobj_bd_SparseCoLoMosek,9))])
disp(' ')

disp(['unstruct: ',num2str(round(fobj_unstruct,9))])