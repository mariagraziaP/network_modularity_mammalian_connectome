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

indices = loadname(fullfile(out_dir, "ModularityMeasures.mat"));


%%
%% GET COMMUNICATION MEASURES
%%

for s=1:S
    
    inv_wei = 1./squeeze(conn_mat(:,:,s));
    m = logical(triu(ones(size(squeeze(conn_mat(:,:,s)),1)),1));

    tmp_coord = squeeze(coords(:,:,s));
    eucl_dist = squareform(pdist(tmp_coord)); 

    % shortest path efficiency
    spe = 1./distance_wei_floyd(inv_wei);
    corr_spe_ed(s) = corr(spe(m), eucl_dist(m), 'type', 'Spearman');
    
    % negative search information
    nsi = -1.*search_information(squeeze(conn_mat(:,:,s)), inv_wei, 0);
    corr_nsi_ed(s) = corr(nsi(m), eucl_dist(m), 'type', 'Spearman');

    % communicability
    cmy = communicability_wei(squeeze(conn_mat(:,:,s)));
    corr_cmy_ed(s) = corr(cmy(m), eucl_dist(m), 'type', 'Spearman');

end

%% correlation with brain volume

% shortest path efficiency

comm_measures(1).name = 'SPE';
comm_measures(1).val = corr_spe_ed;

[rho, pval] = corr(corr_spe_ed', log10_brain_volume, 'Type', 'Spearman');
comm_measures(1).corr_log10BrVol = [rho, pval];

[rho, pval] = corr(corr_spe_ed', indices.brainSize(4).val, 'Type', 'Spearman');
comm_measures(1).corr_log10GreyMat = [rho, pval];

[rho, pval] = corr(corr_spe_ed', indices.brainSize(6).val, 'Type', 'Spearman');
comm_measures(1).corr_log10WhiteMat = [rho, pval];

% negative search information

comm_measures(2).name = 'NSI';
comm_measures(2).val = corr_nsi_ed;

[rho, pval] = corr(corr_nsi_ed', log10_brain_volume, 'Type', 'Spearman');
comm_measures(2).corr_log10BrVol = [rho, pval];

[rho, pval] = corr(corr_nsi_ed', indices.brainSize(4).val, 'Type', 'Spearman');
comm_measures(2).corr_log10GreyMat = [rho, pval];

[rho, pval] = corr(corr_nsi_ed', indices.brainSize(6).val, 'Type', 'Spearman');
comm_measures(2).corr_log10WhiteMat = [rho, pval];

% communicability

comm_measures(3).name = 'CMY';
comm_measures(3).val = corr_cmy_ed;

[rho, pval] = corr(corr_cmy_ed', log10_brain_volume, 'Type', 'Spearman');
comm_measures(3).corr_log10BrVol = [rho, pval];

[rho, pval] = corr(corr_cmy_ed', indices.brainSize(4).val, 'Type', 'Spearman');
comm_measures(3).corr_log10GreyMat = [rho, pval];

[rho, pval] = corr(corr_cmy_ed', indices.brainSize(6).val, 'Type', 'Spearman');
comm_measures(3).corr_log10WhiteMat = [rho, pval];

%% save

save(fullfile(out_dir, 'CommunicationMeasures.mat'), 'measures')

%% 
%% GET COMMUNICATION MEASURES ON NULL NETWORKS (LATTICE AND RANDOM NETS)
%%

numiter = 100;

for s=1:size(conn_mat,3)
    
    parfor i=1:numiter
        
        disp([num2str(s) ' ' num2str(i)])

        tmpnull_rand = makerandCIJ_und(N, nnz(conn_mat(:,:,s))/2);
        tmpnull_lat = makelatticeCIJ(N, nnz(conn_mat(:,:,s))/2);

        inv_wei_rand = 1./tmpnull_rand;
        inv_wei_lat = 1./tmpnull_lat;

        try
            sp = 1./distance_wei_floyd(inv_wei_rand) ;
            avsp_rand(s,i) = mean(sp(m));
            sp = 1./distance_wei_floyd(inv_wei_lat) ;
            avsp_lat(s,i) = mean(sp(m));
        end

        try
            co = communicability_wei(tmpnull_rand) ;
            avco_rand(s,i) = mean(co(m));
            co = communicability_wei(tmpnull_lat) ;
            avco_lat(s,i) = mean(co(m));
        end

        try
            SI = search_information(tmpnull_rand, inv_wei_rand);
            avsi_rand(s,i) = mean(-SI(m));
            SI = search_information(tmpnull_lat, inv_wei_lat);
            avsi_lat(s,i) = mean(-SI(m));
        end

    
    end
end

norm_sp = (avsp'-mean(avsp_lat,2))./(mean(avsp_rand,2)-mean(avsp_lat,2));
norm_co = (avco'-mean(avco_lat,2))./(mean(avco_rand,2)-mean(avco_lat,2));

save(fullfile('...\output', 'NullRandLat'),...
    'avsp_rand', 'avsp_lat', 'avco_rand', 'avco_lat', 'avsi_rand', 'avsi_lat',...
    'avmt_rand', 'avmt_lat', 'norm_sp', 'norm_co')
