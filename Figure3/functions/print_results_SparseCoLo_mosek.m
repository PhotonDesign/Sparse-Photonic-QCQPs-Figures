global Time
disp(['SparseCoLo + Mosek time INIT: ',num2str(Time.SparseCoLO_init),' s']);
disp(['SparseCoLo + Mosek time SOLVE: ',num2str(Time.SparseCoLO),' s']);
disp(['SparseCoLo + Mosek time TOTAL: ',num2str(time_SparseCoLo),' s']);
disp(' ')

disp(['SparseCoLo + Mosek bound: ',num2str(round(fobj_bd_SparseCoLoMosek,9))])
disp(' ')