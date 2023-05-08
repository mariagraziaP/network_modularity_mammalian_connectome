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


%% get MRCC output

load(fullfile(out_dir, 'MRCCmammals_results.mat'))

G = size(all_comm,3);
gamma_range = gamma_range(1:G);

%% get permutations

nper = 1000;
for n=1:nper
    perm_set(:,n) = randperm(200);
end

NullMeasures.permset = perm_set;


%% module density

for n=1:nper
    for s=1:S
        [~, ~, ~, ~, intra_norm(s,n)] = get_InterIntraModulesDensity(...
            squeeze(conn_mat(:,:,s)), squeeze(cons(perm_set(:,n),s)));
    end
    [rho_brvol(n), pval_brvol(n)] = corr(intra_conn_norm(:,n), log10_brain_volume, 'Type', 'Spearman');
    [rho_grm(n), pval_grm(n)] = corr(intra_conn_norm(:,n), table2array(newsheet(:,10)), 'Type', 'Spearman');
    [rho_wm(n), pval_wm(n)] = corr(intra_conn_norm(:,n), table2array(newsheet(:,11)), 'Type', 'Spearman');
end

NullMeasures.ModularityIndices(1).name = 'ModulesDensity';
NullMeasures.ModularityIndices(1).val = intra_conn_norm;

NullMeasures.ModularityIndices(1).corr_log10BrVol = [rho_brvol; pval_brvol];
NullMeasures.ModularityIndices(1).corr_log10GreyMat = [rho_grm; pval_grm];
NullMeasures.ModularityIndices(1).corr_log10WhiteMat = [rho_wm; pval_wm];

clear rho_brvol rho_grm rho_wm pval_brvol pval_grm pval_wm


%% CC,ED

for n=1:nper
    for i=1:size(agree_mat,3)
        [row, col, val] = get_triu_CC(agree_mat(perm_set(:,n),perm_set(:,n),i));
        for j=1:size(row,1)
            ED(j) = EuclideanDistance([coords(row(j),1,i) coords(col(j),1,i)],...
                [coords(row(j),2,i) coords(col(j),2,i)],...
                [coords(row(j),3,i) coords(col(j),3,i)]);
        end
        [rho_cced(i,n), pval_cced(i,n)] = corr(val, ED', 'type','Spearman');
        clear ED 
    end
    [rho_brvol(n), pval_brvol(n)] = corr(rho_cced(:,n), log10_brain_volume, 'Type', 'Spearman');
    [rho_grm(n), pval_grm(n)] = corr(rho_cced(:,n), table2array(newsheet(:,10)), 'Type', 'Spearman');
    [rho_wm(n), pval_wm(n)] = corr(rho_cced(:,n), table2array(newsheet(:,11)), 'Type', 'Spearman');
end

NullMeasures.ModularityIndices(2).name = 'CCED';
NullMeasures.ModularityIndices(2).val = rho_cced;

NullMeasures.ModularityIndices(2).corr_log10BrVol = [rho_brvol; pval_brvol];
NullMeasures.ModularityIndices(2).corr_log10GreyMat = [rho_grm; pval_grm];
NullMeasures.ModularityIndices(2).corr_log10WhiteMat = [rho_wm; pval_wm];

clear rho_brvol rho_grm rho_wm pval_brvol pval_grm pval_wm


%% inter-hem modules

hemind = [ones(100,1); 2*ones(100,1)];

for n=1:nper
    for s=1:S
        alldens = get_IntraInterHemDensity_fromAgreement(squeeze(agree_mat(perm_set(:,n),perm_set(:,n),s)), hemind);
        ihd(s,n) = squeeze(alldens(1,2));
    end
    [rho_brvol(n), pval_brvol(n)] = corr(ihd(:,n), log10_brain_volume, 'Type', 'Spearman');
    [rho_grm(n), pval_grm(n)] = corr(ihd(:,n), table2array(newsheet(:,10)), 'Type', 'Spearman');
    [rho_wm(n), pval_wm(n)] = corr(ihd(:,n), table2array(newsheet(:,11)), 'Type', 'Spearman');
end

NullMeasures.ModularityIndices(3).name = 'InterHemMod';
NullMeasures.ModularityIndices(3).val = ihd;

NullMeasures.ModularityIndices(3).corr_log10BrVol = [rho_brvol; pval_brvol];
NullMeasures.ModularityIndices(3).corr_log10GreyMat = [rho_grm; pval_grm];
NullMeasures.ModularityIndices(3).corr_log10WhiteMat = [rho_wm; pval_wm];

clear rho_brvol rho_grm rho_wm pval_brvol pval_grm pval_wm

%% CC,cost

for i=1:size(agree_mat,3)
    cost = squeeze(conn_mat(:,:,i)).*squeeze(tract_length(:,:,i));
    cost(isnan(cost)) = 0;
    [~, ~, cost_triu] = get_triu_CC(cost);

    for n=1:nper
        [~, ~, cc_triu] = get_triu_CC(agree_mat(perm_set(:,n),perm_set(:,n),i));
        [rho_ccc(i,n), pval_ccc(i,n)] = corr(cost_triu, cc_triu, 'type','Spearman');
    end
end

for n=1:nper
    [rho_brvol(n), pval_brvol(n)] = corr(rho_ccc(:,n), log10_brain_volume, 'Type', 'Spearman');
    [rho_grm(n), pval_grm(n)] = corr(rho_ccc(:,n), table2array(newsheet(:,10)), 'Type', 'Spearman');
    [rho_wm(n), pval_wm(n)] = corr(rho_ccc(:,n), table2array(newsheet(:,11)), 'Type', 'Spearman');
end

NullMeasures.ModularityIndices(4).name = 'CCcost';
NullMeasures.ModularityIndices(4).val = rho_ccc;

NullMeasures.ModularityIndices(4).corr_log10BrVol = [rho_brvol; pval_brvol];
NullMeasures.ModularityIndices(4).corr_log10GreyMat = [rho_grm; pval_grm];
NullMeasures.ModularityIndices(4).corr_log10WhiteMat = [rho_wm; pval_wm];

clear rho_brvol rho_grm rho_wm pval_brvol pval_grm pval_wm

%%  long distance

for n=1:nper
    for s=1:S
        inter_mod_LD(s,n)  = longDistanceConnections_modules(...
            conn_mat(:,:,s), coords(:,:,s), cons(perm_set(:,n),s));
    end
    [rho_brvol(n), pval_brvol(n)] = corr(inter_mod_LD(:,n), log10_brain_volume, 'Type', 'Spearman');
    [rho_grm(n), pval_grm(n)] = corr(inter_mod_LD(:,n), table2array(newsheet(:,10)), 'Type', 'Spearman');
    [rho_wm(n), pval_wm(n)] = corr(inter_mod_LD(:,n), table2array(newsheet(:,11)), 'Type', 'Spearman');
end

NullMeasures.ModularityIndices(5).name = 'longDistConnMod';
NullMeasures.ModularityIndices(5).val = inter_mod_LD;

NullMeasures.ModularityIndices(5).corr_log10BrVol = [rho_brvol; pval_brvol];
NullMeasures.ModularityIndices(5).corr_log10GreyMat = [rho_grm; pval_grm];
NullMeasures.ModularityIndices(5).corr_log10WhiteMat = [rho_wm; pval_wm];

clear rho_brvol rho_grm rho_wm pval_brvol pval_grm pval_wm


%% weight

load('D:\work\mammalian_connectome\output\ModularityMeasures_Null_LabPerm.mat')

for n=275:nper
    disp(n)
    parfor s=1:S
        [~, ~, cc_triu] = get_triu_CC(agree_mat(NullMeasures.permset(:,n),NullMeasures.permset(:,n),s));
        [~, ~, conn_triu] = get_triu_CC(conn_mat(:,:,s));
        [rhoccw(s,n), pvalccw] = corr(cc_triu, conn_triu, 'type', 'Spearman');
    end
    [rho_brvol(n), pval_brvol(n)] = corr(rhoccw(:,n), log10_brain_volume, 'Type', 'Spearman');
    [rho_grm(n), pval_grm(n)] = corr(rhoccw(:,n), table2array(newsheet(:,10)), 'Type', 'Spearman');
    [rho_wm(n), pval_wm(n)] = corr(rhoccw(:,n), table2array(newsheet(:,11)), 'Type', 'Spearman');
end

NullMeasures.ModularityIndices(6).name = 'CCweight';
NullMeasures.ModularityIndices(6).val = rhoccw;

NullMeasures.ModularityIndices(6).corr_log10BrVol = [rho_brvol; pval_brvol];
NullMeasures.ModularityIndices(6).corr_log10GreyMat = [rho_grm; pval_grm];
NullMeasures.ModularityIndices(6).corr_log10WhiteMat = [rho_wm; pval_wm];

clear rho_brvol rho_grm rho_wm pval_brvol pval_grm pval_wm


%% save

save(fullfile(out_dir, 'ModularityMeasures_Null_LabPerm.mat'), 'NullMeasures')

