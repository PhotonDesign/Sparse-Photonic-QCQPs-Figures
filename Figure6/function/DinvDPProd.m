% just for fields, no deriv's (for now)
function[DDP] = DinvDPProd(lp, lpp1, pol, sz)
DI = DInvMatrix(lp, pol, sz);
D = DMatrix(lpp1, pol, sz);
P = PMatrix(lpp1, sz);
DDP = Times3(DI,D,P);
end