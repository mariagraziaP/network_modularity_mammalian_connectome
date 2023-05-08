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

nreps = 2500;

%% get null communication measures

m = logical(triu(ones(size(conn_mat,1)),1));

for s=1:S

    tmp_coord = squeeze(coords(:,:,s));
    eucl_dist = squareform(pdist(tmp_coord));

    for n = 1:nreps
        
        % null
        [~,rewr] = geombinsurr_partial(conn_mat(:,:,s), tract_length(:,:,s),...
            1, 20, 'quantiles') ;

        inv_wei = 1./squeeze(rewr);

        % shortest path efficiency
        spe = 1./distance_wei_floyd(inv_wei);
        corr_spe_ed(s,n) = corr(spe(m), eucl_dist(m), 'type', 'Spearman');

        % negative search information
        nsi = -1.*search_information(rewr, inv_wei, 0);
        corr_nsi_ed(s,n) = corr(nsi(m), eucl_dist(m), 'type', 'Spearman');

        % communicability
        cmy = communicability_wei(rewr);
        corr_cmy_ed(s,n) = corr(cmy(m), eucl_dist(m), 'type', 'Spearman');

    
    end
end

%% correlation with brain volume

for n=1:nreps

    % shortest path efficiency

    [rho_spebv(n), pval_spebv(n)] = corr(corr_spe_ed(:,n), log10_brain_volume, 'Type', 'Spearman');
    [rho_spegm(n), pval_spegm(n)] = corr(corr_spe_ed(:,n), indices.brainSize(4).val, 'Type', 'Spearman');
    [rho_spewm(n), pval_spewm(n)] = corr(corr_spe_ed(:,n), indices.brainSize(6).val, 'Type', 'Spearman');

    % negative search information

    [rho_nsibv(n), pval_nsibv(n)] = corr(corr_nsi_ed(:,n), log10_brain_volume, 'Type', 'Spearman');
    [rho_nsigm(n), pval_nsigm(n)] = corr(corr_nsi_ed(:,n), indices.brainSize(4).val, 'Type', 'Spearman');
    [rho_nsiwm(n), pval_nsiwm(n)] = corr(corr_nsi_ed(:,n), indices.brainSize(6).val, 'Type', 'Spearman');

    % communicability

    [rho_cmybv(n), pval_cmybv(n)] = corr(corr_cmy_ed(:,n), log10_brain_volume, 'Type', 'Spearman');
    [rho_cmygm(n), pval_cmygm(n)] = corr(corr_cmy_ed(:,n), indices.brainSize(4).val, 'Type', 'Spearman');
    [rho_cmywm(n), pval_cmywm(n)] = corr(corr_cmy_ed(:,n), indices.brainSize(6).val, 'Type', 'Spearman');

end

% spe

null_comm_measures(1).name = 'SPE';
null_comm_measures(1).val = corr_spe_ed;

null_comm_measures(1).corr_log10BrVol = [rho_spebv; pval_spebv];
null_comm_measures(1).corr_log10GreyMat = [rho_spegm; pval_spegm];
null_comm_measures(1).corr_log10WhiteMat = [rho_spewm; pval_spewm];

% nsi

null_comm_measures(2).name = 'NSI';
null_comm_measures(2).val = corr_nsi_ed;

null_comm_measures(2).corr_log10BrVol = [rho_nsibv; pval_nsibv];
null_comm_measures(2).corr_log10GreyMat = [rho_nsigm; pval_nsigm];
null_comm_measures(2).corr_log10WhiteMat = [rho_nsiwm; pval_nsiwm];

% cmy

null_comm_measures(3).name = 'CMY';
null_comm_measures(3).val = corr_cmy_ed;

null_comm_measures(3).corr_log10BrVol = [rho_cmybv; pval_cmybv];
null_comm_measures(3).corr_log10GreyMat = [rho_cmygm; pval_cmygm];
null_comm_measures(3).corr_log10WhiteMat = [rho_cmywm; pval_cmywm];


%% save

save(fullfile(out_dir, 'CommunicationMeasures_Null.mat'), 'NullMeasures')
