% pre-requisite: 
% 1. CVX (https://cvxr.com/cvx/)
% 2. Mosek (https://www.mosek.com/)
% 3. SparseCoLO (http://www.opt.c.titech.ac.jp/kojima/SparseCoLO/SparseCoLO.htm)


close all; clc; clear all;
set(0,'DefaultLineLineWidth', 1);
set(0,'DefaultAxesTitleFontWeight','normal');

%% local
cvx_solver mosek
cvx_save_prefs
addpath(genpath('C:\Program Files\Mosek'))
addpath(genpath('functions'))
addpath(genpath('functions_VIE'))

%% cluster
% run /gpfs/loomis/project/owen_miller/sg888/Software/cvx/cvx_startup.m % cvx
% addpath(genpath('/gpfs/loomis/project/owen_miller/sg888/Software/SparseCoLO112/')) % sparsecolo
% addpath(genpath('/gpfs/loomis/project/owen_miller/sg888/Software/sedumi-master/')) % sedumi
% addpath(genpath('/gpfs/loomis/project/owen_miller/sg888/Software/mosek/9.3/toolbox/r2015a/')) % mosek
% addpath(genpath('functions'))
% addpath(genpath('functions_VIE'))

%%
data_folder = 'data';
folder_name = ['sweep_metasurf_focus_',datestr(now,'dd-mmm-yyyy HH.MM.SS'),'/']; % for saving data
mkdir(data_folder, [folder_name,'/'])

%% sweep parameters
lambda = .7; % wavelength (c=1);

chi =1;% 1.25; % material
h = 0.03; % discretization step size
Lx = 0.4;0.875;

thickness = 0.15; %pml thickness 0.3;%
beta = 10; %pml strength 7;%

Ly_list = [0.5100    0.6600    0.8100    0.9600]; [0.6,1,1.5,3];

%calculation mode
const_NA_or_d = 'NA'; %'NA'; %'d'

d_list = 0; % only used when: const_NA_or_d = 'd'
NA_list = 0.9;[0.6,0.75,0.9,0.99]; % only used when: const_NA_or_d = 'NA'

if (strcmp(const_NA_or_d,'NA') && length(d_list)>1 ) ||...
        (strcmp(const_NA_or_d,'d') && length(NA_list)>1)
disp('ERROR: Only d or NA can be swept')
return;
end

% angles to sweep
Ntheta = 5;
theta_range = pi/3;

%also do VIE calculation for comparison
do_VIE = 1;

%% run sweep
iter_total = length(Ly_list)*length(d_list)*length(NA_list);
iter_count = 1;

for ii_Ly = 1:length(Ly_list)
for ii_d = 1:length(d_list)
for ii_NA = 1:length(NA_list)
    
disp(' ');disp('-------------------------------------------------------------------')
disp(['Ly: ',num2str(ii_Ly),' of ',num2str(length(Ly_list)),...
    ',    d: ',num2str(ii_d),' of ',num2str(length(d_list)),...
    ',    NA: ',num2str(ii_NA),' of ',num2str(length(NA_list)),...
    '    (total: ',num2str(iter_count),' of ',num2str(iter_total),')'])
disp('-------------------------------------------------------------------');
    
Ly = Ly_list(ii_Ly);
d = d_list(ii_d);
NA = NA_list(ii_NA);

file_name = sprintf('Ly=%g__NA=%g__d=%g',Ly,NA,d);

sweep_2D_metasurf; % run calculation

clear S b c D AA BB % these take up a lot of space, no need to save
f = fullfile(data_folder,folder_name,[file_name,'.mat']);
save(f); % save data

iter_count = iter_count+1;
end
end
end

disp(' ');disp('DONE!');