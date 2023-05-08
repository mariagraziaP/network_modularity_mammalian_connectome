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
    sprintf('leng_mat_repani_stack_thr%s_.mat', thresholds{curr_thresh})));
tract_length = data;
clear newsheet

load(fullfile(data_dir,...
    sprintf('con_mat_gn_repani_stack_thr%s_.mat', thresholds{curr_thresh})));
conn_mat = data;
clear data

N = size(conn_mat, 1);
S = size(conn_mat, 3);

brain_volume = table2array(newsheet(:,13));
log10_brain_volume = table2array(newsheet(:,12));

load(fullfile(out_dir, 'MRCCmammals_results.mat'))


%%
%% GET CORRELATIONS BETWEEN COCLASSIFICATION AND COMMUNICATION
%%

m = logical(triu(ones(size(conn_mat,1)),1));

for s=1:size(conn_mat,3)
    
    cocl = squeeze(agree_mat(:,:,s));
    inv_wei = 1./conn_mat(:,:,s);
    
    % coclassification vs shortest path efficiency
    try
        sp = 1./distance_wei_floyd(inv_wei) ;
        [corrsp(s), pvalsp(s)] = corr(sp(m), cocl(m),'type','s') ;
        avsp(s) = mean(sp(m));
    end
    
    % coclassification vs communicability
    try
        co = communicability_wei(conn_mat(:,:,s)) ;
        [corrco(s), pvalco(s)] = corr(co(m), cocl(m),'type','s') ;
        avco(s) = mean(co(m));
    end
    
    % coclassification vs negative search information 
    try
        SI = search_information(conn_mat(:,:,s), inv_wei);
        [corrsi(s), pvalsi(s)] = corr(-SI(m), cocl(m),'type','s') ;
        avsi(s) = mean(-SI(m));
    end

end

%% save

save(fullfile('...\output', 'modularity_vs_communication'),...
    'corrsp', 'pvalsp', 'corrco', 'pvalco', 'corrsi', 'pvalsi')