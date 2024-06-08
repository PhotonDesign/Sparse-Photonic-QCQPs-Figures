% Use optimization toolbox function fmincon to find largest transmission of
% an arbitrary initial layered material

function[mHist,xHist,fc,x,fval] = test3c()

% Initial indice set as random numbers
% const_epsr = [2.4025 13.3225 2.4025 13.3225];
% Number of layers
% num_of_layers = length(const_epsr);
% Set layerThickness
% layerThicknesses = [0.05,0.06,0.05,0.06];
num_of_layers = 80;
const_epsr = zeros(1,num_of_layers);
const_epsr(1:2:num_of_layers) = (2.3+0.03*1i)^2;
const_epsr(2:2:num_of_layers) = 1;
layerThicknesses = 0.1*rand(1,num_of_layers);
pol = {'p'};

mHist = [];
xHist = zeros(0,num_of_layers);
fc = [];
% history.x = [];
% history.R = [];

% Wavelength
lambda = 1;
w = 2 * pi ./ lambda;

% Incidence Angle
theta = 0;

target_R = 1*exp(1i*0*pi);

options = optimoptions('fmincon', 'SpecifyObjectiveGradient', true, 'OutputFcn',@outfun, ...
    'Display','iter-detailed','PlotFcn','optimplotfval','CheckGradient',false,'Algorithm','sqp');
fx = @(layerThicknesses) fom(const_epsr, layerThicknesses, w, theta, pol, target_R);
[x,fval] = fmincon(fx,layerThicknesses,[],[],[],[],zeros(num_of_layers,1),(0.1+1e-5)*ones(num_of_layers,1),[],options);
% maxIter = 500;
% [x,fval] = grad_descent(layerThicknesses, const_epsr, w, theta, pol, target_R, maxIter);

lm = createLayerMaterials(const_epsr, x);
[R,~,r] = multilayer_film_derivs_thickness(lm, x, w, theta, pol);
figure; plot(theta, abs(R{1}),'o-')
hold on; plot(theta, abs(target_R),'o-')

    function[x,fval] = grad_descent(x, eps_array, w, theta, pol, target_R, maxIter)
%     thickChangePerIter = 0.005;
    
    tStart = 0.05;
    fc = 0;
    for i = 1:maxIter
        xHist(end+1,:) = x;
        [f, dfdthick, Ri] = fom(eps_array, x, w, theta, pol, target_R);
        mHist(end+1) = Ri{1};
        fc(end+1) = fc(end)+2;
        
        needStep = 1;
        t = tStart;
        beta = 0.5;
        dfdx = dfdthick/sum(abs(dfdthick));
        while (needStep)
            xnew = x - t * dfdx;
            xnew(xnew<0)=0;
            [fnew,~,Rinew] = fom(eps_array, xnew, w, theta, pol, target_R);
            xHist(end+1,:) = xnew;
            mHist(end+1) = Rinew{1};
            fc(end+1) = fc(end)+2;
            if (fnew > (f - t/2 * norm(dfdx)^2))
                t = beta * t;
            else
                needStep = 0;
                x = xnew;
                if (t<0.5*tStart)
                    tStart = 0.5*tStart;
                end
            end
        end
        % without any line search
        % x = x - thickChangePerIter * dfdthick/sum(abs(dfdthick));
        [i f t]
        if (t<1e-4)
            break;
        end
    end
    fval = mHist(end);
    end

    function stop = outfun(x,optimValues,state)
    stop = false;
    switch state
        case 'iter'
%             history.x = [history.x; x];
            layerMaterials = createLayerMaterials(const_epsr, x);
            Ri = multilayer_film_derivs_thickness(layerMaterials, x, w, theta, pol);
%             history.R = [history.R; R{1}];
            mHist(end+1) = Ri{1};
            xHist(end+1,:) = x;
            fc(end+1) = optimValues.funccount;
        otherwise
    end
    end

% dTdn over different w, pol requires additional cycles
% ONE w, MULTIPLE thetas
    function [f, dfdthick, R] = fom(eps_array, layerThicknesses, w, theta, pol, target_R)
    layerMaterials = createLayerMaterials(eps_array, layerThicknesses);
    [R,dR,r,dr,T,dT,t,dt] = multilayer_film_derivs_thickness(layerMaterials, layerThicknesses, w, theta, pol);

    %     % R - tar_R
    %     deltaR = reshape(cell2mat(R),1,length(theta))-target_R;
    %     f = sum(deltaR.^2, 'all');
    %     dR = reshape(cell2mat(dR),length(theta),length(eps_array));
    %     dfdthick = 2*deltaR*dR;
    %
    % r - tar_r
    deltar = reshape(cell2mat(r),1,length(theta))-target_R;
    f = sum(abs(deltar).^2, 'all');
    dr = reshape(cell2mat(dr),length(theta),length(eps_array));
    dfdthick = 2*real(conj(deltar)*dr);

%     mHist(end+1) = R{1};
%     xHist(end+1,:) = layerThicknesses;
    end
end
% function [f, dfdthick] = phase_constraint(eps_array, layerThicknesses, w, theta, pol, target_Phase)
%     layerMaterials = cell(1,length(eps_array)+2);
%     layerMaterials{1} = @Vacuum;
%     layerMaterials{end} = @Vacuum;
%     for i=1:length(layerThicknesses)
%         layerMaterials{i+1} = @(w)eps_array(i);
%     end
%     [R,dR,r,dr,T,dT,t,dt] = multilayer_film_derivs_thickness(layerMaterials, layerThicknesses, w, theta, pol);
%     % R - tar_R
% %     deltaR = reshape(cell2mat(R),1,length(theta))-target_R;
%     % sum {R- tar_R}^2
% %     f = sum(deltaR.^2, 'all');
%
%     f = angle()
%     % dRdthick
%     dR = reshape(cell2mat(dR),length(theta),length(eps_array));
%     % sum{2 * (R-tar_R) *dR}
%     dfdthick = 2*deltaR*dR;
%
%     x0 = 0;
%     x = linspace(x0,x0+10,200);
%
%     global mHist;
%     mHist(end+1) = R{1};
%     global xHist;
%     xHist(end+1,:) = layerThicknesses;
% end



% % dTdn over different theta, w, pol requires additional cycles
% function [T, dT] = fom(eps_array, layerThicknesses, w, theta, pol)
%     layerMaterials = cell(1,length(eps_array)+2);
%     layerMaterials{1} = @Vacuum;
%     layerMaterials{end} = @Vacuum;
%     for i=1:length(layerThicknesses)
%         layerMaterials{i+1} = @(w)eps_array(i);
%     end
%     [R,dR,r,dr,T,dT,t,dt] = multilayer_film_derivs_thickness(layerMaterials, layerThicknesses, w, theta, pol);
%     T = 1.0 .* cell2mat(T);
%     dT = 1.0 .* reshape(cell2mat(dT),length(eps_array),1);
% end
