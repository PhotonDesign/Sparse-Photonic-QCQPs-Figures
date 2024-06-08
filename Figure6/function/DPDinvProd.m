% dFlag = 1 if want deriv's
function[DPD, dDPDdn] = DPDinvProd(lp, pol, sz, dFlag)
D = DMatrix(lp, pol, sz);
DI = DInvMatrix(lp, pol, sz);
P = PMatrix(lp, sz);
DPD = Times3(D,P,DI);

if (dFlag)
    dDdn = dDMatrix(lp, pol, sz);
    dDIdn = dDInvMatrix(lp, pol, sz);
    dPdn = dPMatrix(lp, sz);
    dDPDdn = Times3(dDdn,P,DI) + Times3(D,dPdn,DI) + Times3(D,P,dDIdn);
end
end