
function[mHist,xHist,fc,xval,r] = inverse_designMMA(alpha0)

k = 2*pi;
xDesignStart = 0;
L = 4; % size in wavelengths
Nx = 100;
nMat = 2.3+0.03i;
rTar = exp(0.3*pi*1i);
betaStart = 1; % artificial loss parameter to promote binary materials
betaEnd = 10;
beta = betaStart;
maxMMAoutiter  = 200; 

[x,dx,Nx,designInd,Ndesign] = set_fdfd_grid(L,Nx,xDesignStart);

if (nargin < 1)
    % % % near-zero start
    alpha0 = 0.1*rand(Ndesign,1); % initial density, set near 0 (start from "scratch")
    
%     % % % random-thickness start
%     w = 0.1*rand(1,80);
%     x0 = xDesignStart;
%     xDesign = x(designInd);
%     alpha0 = zeros(Ndesign,1);
%     for i=1:80
%         xEnd = x0 + w(i);
%         if (mod(i,2)==1) % every other layer is high index
%             alpha0(logical((xDesign>x0).*(xDesign<=xEnd))) = 1;
%         end
%         x0 = xEnd;
%     end
end

lb = zeros(Ndesign,1); % density bounds
ub = ones(Ndesign,1);

warning('off','MATLAB:singularMatrix')
warning('off','MATLAB:nearlySingularMatrix')

m = 1;
n = Ndesign;
epsimin = 0.0000001;
xval    = alpha0;
xold1   = xval;
xold2   = xval;
xmin    = lb;
xmax    = ub;
low     = xmin;
upp     = xmax;
c       = zeros(m,1);
d       = zeros(m,1);
a0      = 0;
a       = zeros(m,1);
outeriter = 0;
kkttol  = 0;


mHist = [];
xHist = zeros(0,Ndesign);
fc = [];

figure(2); plot(x(designInd),xval);
[f0val,df0dx,fval,dfdx,r] = sim_dir_adj(xval,x,Nx,dx,designInd,nMat,k,rTar,beta);
[outeriter f0val]

fc(1) = 1;
mHist(end+1) = abs(r)^2;
xHist(end+1,:) = xval;

%%%% The iterations start
if (nargin < 1)
    kktnorm = kkttol+10;
    outit = 0;
    while kktnorm > kkttol && outit < maxMMAoutiter
        outit   = outit+1;
        outeriter = outeriter+1;

        [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp] = ...
            mmasub(m,n,outeriter,xval,xmin,xmax,xold1,xold2, ...
            f0val,df0dx,fval,dfdx,low,upp,a0,a,c,d);

        xold2 = xold1;
        xold1 = xval;
        xval  = xmma;

        figure(2); plot(x(designInd),xval);
        [f0val,df0dx,fval,dfdx,r] = sim_dir_adj(xval,x,Nx,dx,designInd,nMat,k,rTar,beta);
        if(mod(outit,10)==0)
            beta = betaStart + 2*outit/maxMMAoutiter * (betaEnd-betaStart);
            if (beta>betaEnd), beta = betaEnd; end
        end

        fc(end+1) = fc(end)+1;
        mHist(end+1) = abs(r)^2;
        xHist(end+1,:) = xval;

        [residu,kktnorm,residumax] = ...
            kktcheck(m,n,xmma,ymma,zmma,lam,xsi,eta,mu,zet,s, ...
            xmin,xmax,df0dx,fval,dfdx,a0,a,c,d);
        [outeriter f0val abs(r)^2]
        %
    end
end

end