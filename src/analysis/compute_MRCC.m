%%
clear
clc

%% set path

addpath('...\fcn')
addpath(genpath('...\ext'))

%% set directories

data_dir = '...\mammals_connectome';
out_dir = '...\output';

%% load data and set parameters

thresholds = {'0' '0.05' '0.1' '0.15'};
curr_thresh = 1;

load(fullfile(data_dir,...
    sprintf('cents_repani_stack_thr%s_.mat', thresholds{curr_thresh})));
coords = data;
clear data newsheet

load(fullfile(data_dir,...
    sprintf('con_mat_gn_repani_stack_thr%s_.mat', thresholds{curr_thresh})));

N = size(data, 1);
S = size(data, 3);

%% compute MRCC new function get_MRCC_mammals - normilize conn mat

gamma_sample = 500;
in_gamma_low = 0.01;
in_gamma_high = 10;

for s=1:S
    mat_norm = squeeze(data(:,:,s));
    mat_norm = mat_norm./(max(mat_norm(:)));
    mat_norm = mat_norm.*~eye(size(mat_norm,1));
    data_norm(:,:,s) = mat_norm;
end

[cons, agree_mat, anull, A, all_comm, gamma_range] = get_MRCC_mammals(...
    data_norm, gamma_sample, gamma_sample, 100, in_gamma_low, in_gamma_high, 0);

save(fullfile(out_dir, 'MRCCmammals_results.mat'),...
    'cons', 'agree_mat', 'anull', 'A', 'all_comm', 'gamma_range')