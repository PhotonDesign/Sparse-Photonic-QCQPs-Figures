function [Pref,Eref,phase] = getPrefFromChiVac(chi,L,M)
    c = 1;
    wl = 1;
    f = c/wl;
    w = 2*pi*f;
    % Units normalized to 1
    zMat = linspace(0,L,M)';
    % Pick desired length and discretization
    
    %Break into individual blocks
    %Number of breaks = number of interior regions of vacuum
    blockInds = find(chi > 1e-4);
    breaks = find(ischange(blockInds,'linear','threshold',2));
    blocks = cell(length(breaks)+1,1);
    z0 = 1;
    for i = 1:length(breaks)
        blocks{i} = blockInds(z0:(breaks(i)-1));
        z0 = breaks(i);
    end
    blocks{end} = blockInds(breaks(end):length(blockInds));
    
    %Solve for each p-current individually
    %Then sum integrals to get reflected field
    EMat = get_E(zMat,2*pi);
    EMats = cell(length(blocks),1);
    chiMats = cell(length(blocks),1);
    GMats = cell(length(blocks),1);
    phiMats = cell(length(blocks),1);
    N = 1000;
    zExt = linspace(-1,0,N);
    Eref_c = zeros(N,1);
    GMatExts = cell(length(blocks),1);
    for i = 1:length(blocks)
        EMats{i} = EMat(blocks{i});
        chiMats{i} = diag(1./chi(blocks{i}));
        GMats{i} = zeros(length(blocks{i}));
        for a = 1:length(blocks{i})
            for b = 1:length(blocks{i})
                GMats{i}(a,b) = get_G(zMat(blocks{i}(a)),zMat(blocks{i}(b)),2*pi,L/M);
            end
        end
        phiMats{i} = linsolve(GMats{i}-chiMats{i},-EMats{i});
        GMatExts{i} = zeros(N,length(blocks{i})); %Higher resolution grid to account for phase
        for a = 1:N
            for b = 1:length(blocks{i})
                GMatExts{i}(a,b) = get_G(zExt(a),zMat(blocks{i}(b)),2*pi,L/M);
            end
        end
        Eref_c = Eref_c + GMatExts{i}*phiMats{i};
    end
    
    Pref = abs(Eref_c(1))^2;
    Eref = real(Eref_c);
    [~,maxInd] = max(Eref);
    phase = -zExt(maxInd)*2*pi;
    if phase > pi
        phase = phase-2*pi;
    end
end