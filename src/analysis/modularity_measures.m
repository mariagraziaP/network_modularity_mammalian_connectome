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

%% get all measures and store them in a struct

%% data

measures.brainSize(1).name = 'Brain Volume';
measures.brainSize(1).val = brain_volume;

measures.brainSize(2).name = 'log10 Brain Volume';
measures.brainSize(2).val = log10_brain_volume;

measures.brainSize(3).name = 'Grey Matter';
measures.brainSize(3).val = table2array(newsheet(:,8));

measures.brainSize(4).name = 'log10 Grey Matter';
measures.brainSize(4).val = table2array(newsheet(:,10));

measures.brainSize(5).name = 'White Matter';
measures.brainSize(5).val = table2array(newsheet(:,9));

measures.brainSize(6).name = 'log10 White Matter';
measures.brainSize(6).val = table2array(newsheet(:,11));

%%
%% GET MEASURES ON CONSENSUS/AGREEMENT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% module density

for s=1:S
    [~, ~, ~, ~, intra_norm(s)] = get_InterIntraModulesDensity(...
        squeeze(conn_mat(:,:,s)), squeeze(cons(:,s)));
end

measures.ModularityIndices(1).name = 'ModulesDensity';
measures.ModularityIndices(1).val = intra_norm';

[rho, pval] = corr(intra_norm', log10_brain_volume, 'Type', 'Spearman');
measures.ModularityIndices(1).corr_log10BrVol = [rho, pval];

[rho, pval] = corr(intra_norm', measures.brainSize(4).val, 'Type', 'Spearman');
measures.ModularityIndices(1).corr_log10GreyMat = [rho, pval];

[rho, pval] = corr(intra_norm', measures.brainSize(6).val, 'Type', 'Spearman');
measures.ModularityIndices(1).corr_log10WhiteMat = [rho, pval];


%% CC,ED

for i=1:size(agree_mat,3)
    [row, col, val] = get_triu_CC(agree_mat(:,:,i));
    for j=1:size(row,1)
        ED(j) = EuclideanDistance([coords(row(j),1,i) coords(col(j),1,i)],...
            [coords(row(j),2,i) coords(col(j),2,i)],...
            [coords(row(j),3,i) coords(col(j),3,i)]);
    end
    [rho_cced(i), pval_cced(i)] = corr(val, ED', 'type','Spearman');
    clear ED 
end

measures.ModularityIndices(2).name = 'CCED';
measures.ModularityIndices(2).val = rho_cced';

[rho, pval] = corr(rho_cced', log10_brain_volume, 'Type', 'Spearman');
measures.ModularityIndices(2).corr_log10BrVol = [rho, pval];

[rho, pval] = corr(rho_cced', measures.brainSize(4).val, 'Type', 'Spearman');
measures.ModularityIndices(2).corr_log10GreyMat = [rho, pval];

[rho, pval] = corr(rho_cced', measures.brainSize(6).val, 'Type', 'Spearman');
measures.ModularityIndices(2).corr_log10WhiteMat = [rho, pval];

%% inter-hem modules

hemind = [ones(100,1); 2*ones(100,1)];

for s=1:S
    alldens(:,:,s) = get_IntraInterHemDensity_fromAgreement(squeeze(agree_mat(:,:,s)), hemind);
    ihd(s) = squeeze(alldens(1,2,s));
    ihdn(s) = ihd(s)./(mean([squeeze(alldens(1,1,s)), squeeze(alldens(2,2,s))]));
end

measures.ModularityIndices(3).name = 'InterHemMod';
measures.ModularityIndices(3).val = ihd';

[rho, pval] = corr(ihd', log10_brain_volume, 'Type', 'Spearman');
measures.ModularityIndices(3).corr_log10BrVol = [rho, pval];

[rho, pval] = corr(ihd', measures.brainSize(4).val, 'Type', 'Spearman');
measures.ModularityIndices(3).corr_log10GreyMat = [rho, pval];

[rho, pval] = corr(ihd', measures.brainSize(6).val, 'Type', 'Spearman');
measures.ModularityIndices(3).corr_log10WhiteMat = [rho, pval];

%% CC,cost

for i=1:size(agree_mat,3)
    cost = squeeze(conn_mat(:,:,i)).*squeeze(tract_length(:,:,i));
    cost(isnan(cost)) = 0;
    [~, ~, cost_triu] = get_triu_CC(cost);
    [~, ~, cc_triu] = get_triu_CC(agree_mat(:,:,i));

    [rho_ccc(i), pval_ccc(i)] = corr(cost_triu, cc_triu, 'type','Spearman');
end

measures.ModularityIndices(4).name = 'CCcost';
measures.ModularityIndices(4).val = rho_ccc';

[rho, pval] = corr(rho_ccc', log10_brain_volume, 'Type', 'Spearman');
measures.ModularityIndices(4).corr_log10BrVol = [rho, pval];

[rho, pval] = corr(rho_ccc', measures.brainSize(4).val, 'Type', 'Spearman');
measures.ModularityIndices(4).corr_log10GreyMat = [rho, pval];

[rho, pval] = corr(rho_ccc', measures.brainSize(6).val, 'Type', 'Spearman');
measures.ModularityIndices(4).corr_log10WhiteMat = [rho, pval];


%% long distance

for s=1:S
    inter_mod_LD(s)  = longDistanceConnections_modules(...
        conn_mat(:,:,s), coords(:,:,s), cons(:,s));
end

measures.ModularityIndices(5).name = 'longDistConnMod';
measures.ModularityIndices(5).val = inter_mod_LD';

[rho, pval] = corr(inter_mod_LD', log10_brain_volume, 'Type', 'Spearman');
measures.ModularityIndices(5).corr_log10BrVol = [rho, pval];

[rho, pval] = corr(inter_mod_LD', measures.brainSize(4).val, 'Type', 'Spearman');
measures.ModularityIndices(5).corr_log10GreyMat = [rho, pval];

[rho, pval] = corr(inter_mod_LD', measures.brainSize(6).val, 'Type', 'Spearman');
measures.ModularityIndices(5).corr_log10WhiteMat = [rho, pval];


%% weight

measures = loadname(fullfile(out_dir, 'ModularityMeasures.mat'));

for s=1:S
    [~, ~, cc_triu] = get_triu_CC(agree_mat(:,:,s));
    [~, ~, conn_triu] = get_triu_CC(conn_mat(:,:,s));
    [rhoccw(s), pvalccw(s)] = corr(cc_triu, conn_triu, 'type', 'Spearman');
end

measures.ModularityIndices(6).name = 'CCweight';
measures.ModularityIndices(6).val = rhoccw';

[rho, pval] = corr(rhoccw', log10_brain_volume, 'Type', 'Spearman');
measures.ModularityIndices(6).corr_log10BrVol = [rho, pval];

[rho, pval] = corr(rhoccw', measures.brainSize(4).val, 'Type', 'Spearman');
measures.ModularityIndices(6).corr_log10GreyMat = [rho, pval];

[rho, pval] = corr(rhoccw', measures.brainSize(6).val, 'Type', 'Spearman');
measures.ModularityIndices(6).corr_log10WhiteMat = [rho, pval];

%% modules dimension and number of modules

for s=1:S
    num_mod(s) = length(unique(cons(:,s)));
    labels = unique(cons(:,s)); 
    for j=1:length(labels)
        tmpdim(j) = length(find(cons(:,s)==labels(j)));
    end
    mod_dim(s)  = median(tmpdim);
    clear tmpdim
end

measures.ModularityIndices(7).name = 'NumClust';
measures.ModularityIndices(7).val = mod_dim';

[rho, pval] = corr(mod_dim', log10_brain_volume, 'Type', 'Spearman');
measures.ModularityIndices(7).corr_log10BrVol = [rho, pval];

[rho, pval] = corr(rhoccw', measures.brainSize(4).val, 'Type', 'Spearman');
measures.ModularityIndices(7).corr_log10GreyMat = [rho, pval];

[rho, pval] = corr(rhoccw', measures.brainSize(6).val, 'Type', 'Spearman');
measures.ModularityIndices(7).corr_log10WhiteMat = [rho, pval];


%%
%% GET MEASURES FOR EVERY GAMMA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% module density

for s=1:S
    for g=1:G
        [~, ~, ~, ~, intra_norm(s)] = get_InterIntraModulesDensity(squeeze(...
            conn_mat(:,:,s)), squeeze(all_comm(:,s,g)));
    end
end

for g=1:G
    [rho_brvol(g), pval_brvol(g)] = corr(...
        intra_norm(:,g), log10_brain_volume, 'type','Spearman');

    [rho_grmat(g), pval_grmat(g)] = corr(...
        intra_norm(:,g), measures.brainSize(4).val, 'Type', 'Spearman');

    [rho_whmat(g), pval_whmat(g)] = corr(...
        intra_norm(:,g), measures.brainSize(6).val, 'Type', 'Spearman');

end

measures.ModularityIndices_allG(1).name = 'ModulesDensity';
measures.ModularityIndices_allG(1).val = intra_conn_norm;

measures.ModularityIndices_allG(1).corr_log10BrVol = [rho_brvol; pval_brvol];
measures.ModularityIndices_allG(1).corr_log10GreyMat = [rho_grmat; pval_grmat];
measures.ModularityIndices_allG(1).corr_log10WhiteMat = [rho_whmat; pval_whmat];


%% ED

for s=1:S
    for g=1:G
        ed_ratio(s,g) = get_EDratio(conn_mat(:,:,s), coords(:,:,s), squeeze(all_comm(:,s,g)));
    end
end

for g=1:G
    mask = find(not(isnan(ed_ratio(:,g))));
    if length(mask)>1
        [rho_brvol(g), pval_brvol(g)] = corr(...
            ed_ratio(mask,g), log10_brain_volume(mask), 'type','Spearman');
    
        [rho_grmat(g), pval_grmat(g)] = corr(...
            ed_ratio(mask,g), measures.brainSize(4).val(mask), 'Type', 'Spearman');
    
        [rho_whmat(g), pval_whmat(g)] = corr(...
            ed_ratio(mask,g), measures.brainSize(6).val(mask), 'Type', 'Spearman');
    else
        rho_brvol(g) = NaN; pval_brvol(g) = NaN;
    
        rho_grmat(g) = NaN; pval_grmat(g) = NaN;
    
        rho_whmat(g) = NaN; pval_whmat(g) = NaN;
    end
end

measures.ModularityIndices_allG(2).name = 'CCED';
measures.ModularityIndices_allG(2).val = ed_ratio;

measures.ModularityIndices_allG(2).corr_log10BrVol = [rho_brvol; pval_brvol];
measures.ModularityIndices_allG(2).corr_log10GreyMat = [rho_grmat; pval_grmat];
measures.ModularityIndices_allG(2).corr_log10WhiteMat = [rho_whmat; pval_whmat];


%% inter-hem modules

LH = 1:100;
RH = 101:200;

for s=1:S
    for g=1:G
        tmp_labels = squeeze(all_comm(:,s,g));
        
        count(s,g) = 0;
        for i=1:length(tmp_labels)
            pos = find(squeeze(all_comm(:,s,g))==tmp_labels(i));
            if length(find(ismember(LH,pos)))>0 && length(find(ismember(RH,pos)))>0
                count(s,g) = count(s,g)+1;
            end
        end
        ihd(s,g) = count(s,g)/length(tmp_labels);
    end
end

for g=1:G
    [rho_brvol(g), pval_brvol(g)] = corr(...
        ihd(:,g), log10_brain_volume, 'type','Spearman');

    [rho_grmat(g), pval_grmat(g)] = corr(...
        ihd(:,g), measures.brainSize(4).val, 'Type', 'Spearman');

    [rho_whmat(g), pval_whmat(g)] = corr(...
        ihd(:,g), measures.brainSize(6).val, 'Type', 'Spearman');

end

measures.ModularityIndices_allG(3).name = 'InterHemMod';
measures.ModularityIndices_allG(3).val = ihd;

measures.ModularityIndices_allG(3).corr_log10BrVol = [rho_brvol; pval_brvol];
measures.ModularityIndices_allG(3).corr_log10GreyMat = [rho_grmat; pval_grmat];
measures.ModularityIndices_allG(3).corr_log10WhiteMat = [rho_whmat; pval_whmat];

%% cost

for s=1:S
    cost = squeeze(conn_mat(:,:,s)).*squeeze(tract_length(:,:,s));
    for g=1:G
        cost_ratio(s,g) = get_costRatio(cost, squeeze(all_comm(:,s,g)));
    end
end

for g=1:G
    [rho_brvol(g), pval_brvol(g)] = corr(...
        cost_ratio(:,g), log10_brain_volume, 'type','Spearman');

    [rho_grmat(g), pval_grmat(g)] = corr(...
        cost_ratio(:,g), measures.brainSize(4).val, 'Type', 'Spearman');

    [rho_whmat(g), pval_whmat(g)] = corr(...
        cost_ratio(:,g), measures.brainSize(6).val, 'Type', 'Spearman');
end

measures.ModularityIndices_allG(4).name = 'CCcost';
measures.ModularityIndices_allG(4).val = cost_ratio;

measures.ModularityIndices_allG(4).corr_log10BrVol = [rho_brvol; pval_brvol];
measures.ModularityIndices_allG(4).corr_log10GreyMat = [rho_grmat; pval_grmat];
measures.ModularityIndices_allG(4).corr_log10WhiteMat = [rho_whmat; pval_whmat];


%% long distance

for s=1:S
    for g=1:G
        inter_mod_LD(s,g)  = longDistanceConnections_modules(...
            conn_mat(:,:,s), coords(:,:,s), all_comm(:,s,g));
    end
end

for g=1:G
    [rho_brvol(g), pval_brvol(g)] = corr(...
        inter_mod_LD(:,g), log10_brain_volume, 'type','Spearman');

    [rho_grmat(g), pval_grmat(g)] = corr(...
        inter_mod_LD(:,g), measures.brainSize(4).val, 'Type', 'Spearman');

    [rho_whmat(g), pval_whmat(g)] = corr(...
        inter_mod_LD(:,g), measures.brainSize(6).val, 'Type', 'Spearman');
end

measures.ModularityIndices_allG(5).name = 'longDistConnMod';
measures.ModularityIndices_allG(5).val = inter_mod_LD;

measures.ModularityIndices_allG(5).corr_log10BrVol = [rho_brvol; pval_brvol];
measures.ModularityIndices_allG(5).corr_log10GreyMat = [rho_grmat; pval_grmat];
measures.ModularityIndices_allG(5).corr_log10WhiteMat = [rho_whmat; pval_whmat];


%% weight

for s=1:S
    for g=1:G
        wei_ratio(s,g) = get_weightRatio(conn_mat(:,:,s), squeeze(all_comm(:,s,g)));
    end
end

for g=1:G
    [rho_brvol(g), pval_brvol(g)] = corr(...
        wei_ratio(:,g), log10_brain_volume, 'type','Spearman');

    [rho_grmat(g), pval_grmat(g)] = corr(...
        wei_ratio(:,g), measures.brainSize(4).val, 'Type', 'Spearman');

    [rho_whmat(g), pval_whmat(g)] = corr(...
        wei_ratio(:,g), measures.brainSize(6).val, 'Type', 'Spearman');
end

measures.ModularityIndices_allG(6).name = 'CCweight';
measures.ModularityIndices_allG(6).val = wei_ratio;

measures.ModularityIndices_allG(6).corr_log10BrVol = [rho_brvol; pval_brvol];
measures.ModularityIndices_allG(6).corr_log10GreyMat = [rho_grmat; pval_grmat];
measures.ModularityIndices_allG(6).corr_log10WhiteMat = [rho_whmat; pval_whmat];


%% save measures

save(fullfile(out_dir, 'ModularityMeasures.mat'), 'measures')
