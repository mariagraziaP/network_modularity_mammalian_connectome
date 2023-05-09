clear
clc

%% set path

addpath('...\fcn')
addpath(genpath('...\ext'))

%% set directories

data_dir = '...\mammals_connectome';
out_dir = '...\output';

%% load data and set parameters

N = 200;
S = 201;

order = loadname(fullfile(data_dir, "AnimalOrder.mat"));
animal_name = loadname(fullfile(data_dir, "AnimalName.mat"));
order_index = loadname(fullfile(data_dir, "OrderIndex.mat"));

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

load(fullfile(out_dir, 'MRCCmammals_results.mat'))

mod_ind = loadname(fullfile(out_dir, "ModularityMeasures.mat"));
comm_ind = loadname(fullfile(out_dir, "CommunicationMeasures.mat"));
load(fullfile(out_dir, 'modularity_vs_communication.mat'))
load(fullfile(out_dir, 'NullRandLat'))

colors = turbo(12);
ordcol = [7 11 9 1 6 3 2 5 8 12 10 4];
colors = colors(ordcol,:);

col_points = zeros(S, 3);
for odx=1:length(order)
    pos = find(order_index==odx);
    col_points(pos,:) = repmat(colors(odx,:), length(pos), 1);
end

x_lim = [1 7];
x_tick = 2:6;


%% FIGURE 2: BRAIN VOLUME vs WHITE MATTER / GREY MATTER VOLUME

for odx=1:length(order)
    pos = find(order_index==odx);
    if length(pos)>2
        std_brVol(odx) = std(mod_ind.brainSize(2).val(pos));
    end
    data2plot{odx} = mod_ind.brainSize(2).val(pos);
    col2plot(odx,:,:) = cat(3, colors(odx,:), [0 0 0]);
    mean_brvol(odx) = mean(mod_ind.brainSize(2).val(pos));
end
[mean_brvol_ord, mean_ord] = sort(mean_brvol, 'ascend');

% find correlations slope and intercept of the first two graphs
Nspecies = 201;
nboot = 10000;
for i=1:nboot
    ind = randperm(Nspecies, 181);
    % br vol vs gr mat
    [slope_vg(i), intercept_vg(i)] = get_regress_par(mod_ind.brainSize(2).val(ind), mod_ind.brainSize(4).val(ind));
    % gr mat vs wh mat
    [slope_gw(i), intercept_gw(i)] = get_regress_par(mod_ind.brainSize(4).val(ind), mod_ind.brainSize(6).val(ind));
end
av_slope_vg = mean(slope_vg); std_slope_vg = std(slope_vg);
av_slope_gw = mean(slope_gw); std_slope_gw = std(slope_gw);
av_int_vg = mean(intercept_vg); std_int_vg = std(intercept_vg);
av_int_gw = mean(intercept_gw); std_int_gw = std(intercept_gw);
[rho_vg, pval_vg] = corr(mod_ind.brainSize(2).val, mod_ind.brainSize(4).val,...
    'type','Spearman');
[rho_gw, pval_gw] = corr(mod_ind.brainSize(4).val, mod_ind.brainSize(6).val,...
    'type','Spearman');
[rho_vw, pval_vw] = corr(mod_ind.brainSize(2).val, mod_ind.brainSize(6).val,...
    'type','Spearman');

ff = figure;
pp1 = scatter(mod_ind.brainSize(2).val, mod_ind.brainSize(4).val, 50,...
    'MarkerFaceColor', rgb('LightSlateGrey'), 'MarkerEdgeColor', [0 0 0]);
axis equal
xlabel('log_{10} Brain Volume (mm^3)')
ylabel('log_{10} GrayMatter Volume (mm^3)')
grid on; grid minor; box on
set(gca, 'FontSize', 16, 'XLim', [1 7], 'XTick', 2:6, 'YLim', [1 7], 'Ytick', 2:6,...
    'GridColor', [0 0 0], 'MinorGridColor', [0.5 0.5 0.5])
title('\rho = 0.9979, pval < 10^{-15}', 'FontSize', 15.5)
exportgraphics(ff, 'brVol_grMatVol.png', 'Resolution', 300)

ff = figure;
scatter(mod_ind.brainSize(4).val, mod_ind.brainSize(6).val, 50,...
    'MarkerFaceColor', rgb('LightSlateGrey'), 'MarkerEdgeColor', [0 0 0]);
axis equal
xlabel('log_{10} GreyMatter Volume (mm^3)')
ylabel('log_{10} WhiteMatter Volume (mm^3)')
grid on; grid minor; box on
set(gca, 'FontSize', 16, 'XLim', [0 7], 'XTick', 1:6, 'YLim', [0 7], 'Ytick', 1:6,...
    'GridColor', [0 0 0], 'MinorGridColor', [0.5 0.5 0.5])
title('\rho = 0.9644, pval < 10^{-15}', 'FontSize', 15.5)
exportgraphics(ff, 'grMatVol_whMatVol.png', 'Resolution', 300)

ff = figure;
plot_distributions(data2plot(mean_ord), col2plot(mean_ord,:,:));
xlim([0 12.8])
grid on; grid minor; box on
set(gca, 'XTick', 1:length(order), 'XTickLabel', order(mean_ord),...
    'YTick', 2:6, 'FontSize', 16, ...
    'GridColor', [0 0 0], 'MinorGridColor', [0.5 0.5 0.5])
ylabel('log_{10} Brain Volume (mm^3)', 'FontSize', 16)
set(gcf, 'Position', [30  100.0000  480.0000  500.0000])
exportgraphics(ff, 'brVol_taxonomies.png', 'Resolution', 300)


%% FIGURE 3: figures modularity

% mod density
ff = figure;
plot_scatter_mammals(mod_ind.brainSize(2).val, mod_ind.ModularityIndices(1).val,...
        col_points, x_lim, x_tick, [6 20.5], 7:4:19);
xlabel('log_{10} Brain Volume (mm^3)')
ylabel('intra-module density')
title(sprintf('\\rho = %g; pval < 10^{-11}',...
    round(mod_ind.ModularityIndices(1).corr_log10BrVol(1) ,2)))
set(gca, 'FontSize', 16, 'GridColor', [0 0 0],...
    'MinorGridColor', [0.5 0.5 0.5])
box on
exportgraphics(ff, 'mod_density_brvol.png', 'Resolution', 300)

% weight
ff = figure;
plot_scatter_mammals(mod_ind.brainSize(2).val, mod_ind.ModularityIndices(6).val,...
        col_points, x_lim, x_tick, [0.24 0.57], 0.25:0.1:0.55);
ylabel('\rho (CC,weight)')
xlabel('log_{10} Brain Volume (mm^3)')
title(sprintf('\\rho = %g; pval < 10^{-15}',...
    round(mod_ind.ModularityIndices(6).corr_log10BrVol(1) ,2)))
set(gca, 'FontSize', 16, 'GridColor', [0 0 0],...
    'MinorGridColor', [0.5 0.5 0.5])
box on
exportgraphics(ff, 'cc_weight_brvol.png', 'Resolution', 300)

% cost
ff = figure;
plot_scatter_mammals(mod_ind.brainSize(2).val, mod_ind.ModularityIndices(4).val,...
        col_points, x_lim, x_tick, [0.25 0.55], 0.3:0.05:0.5);
xlabel('log_{10} Brain Volume (mm^3)')
ylabel('\rho (CC,cost)')
title(sprintf('\\rho = %g; pval < 10^{-15}',...
    round(mod_ind.ModularityIndices(4).corr_log10BrVol(1) ,2)))
set(gca, 'FontSize', 16, 'GridColor', [0 0 0],...
    'MinorGridColor', [0.5 0.5 0.5])
box on
exportgraphics(ff, 'cc_cost_brvol.png', 'Resolution', 300)

% ed
ff = figure;
plot_scatter_mammals(mod_ind.brainSize(2).val, mod_ind.ModularityIndices(2).val,...
        col_points, x_lim, x_tick, [-0.57 -0.2], -0.55:0.1:-0.25);
xlabel('log_{10} Brain Volume (mm^3)')
ylabel('\rho (CC,ED)')
title(sprintf('\\rho = %g; pval < 10^{-10}',...
    round(mod_ind.ModularityIndices(2).corr_log10BrVol(1) ,2)))
set(gca, 'FontSize', 16, 'GridColor', [0 0 0],...
    'MinorGridColor', [0.5 0.5 0.5])
box on
exportgraphics(ff, 'cc_ed_brvol.png', 'Resolution', 300)

% long distance
ff = figure;
plot_scatter_mammals(mod_ind.brainSize(2).val, mod_ind.ModularityIndices(5).val,...
        col_points, x_lim, x_tick, [0.07 0.142], 0.08:0.02:0.14);
xlabel('log_{10} Brain Volume (mm^3)')
ylabel({'long-dist betw. modules'})
title(sprintf('\\rho = %g; pval < 10^{-6}',...
    round(mod_ind.ModularityIndices(5).corr_log10BrVol(1) ,2)))
set(gca, 'FontSize', 16, 'GridColor', [0 0 0],...
    'MinorGridColor', [0.5 0.5 0.5])
box on
exportgraphics(ff, 'longdist_brvol.png', 'Resolution', 300)

% inter-hemispheric
ff = figure;
plot_scatter_mammals(mod_ind.brainSize(2).val, mod_ind.ModularityIndices(3).val,...
        col_points, x_lim, x_tick, [0.0045 0.042], 0.01:0.01:0.04);
xlabel('log_{10} Brain Volume (mm^3)')
ylabel({'inter-hemispheric modules'})
title(sprintf('\\rho = %g; pval = %g',...
    round(mod_ind.ModularityIndices(3).corr_log10BrVol(1) ,2),...
    round(mod_ind.ModularityIndices(3).corr_log10BrVol(2) ,3)))
set(gca, 'FontSize', 16, 'GridColor', [0 0 0], 'MinorGridColor', [0.5 0.5 0.5])
box on
exportgraphics(ff, 'interhem_brvol.png', 'Resolution', 300)


%% FIGURE 4: communication indices

% shortest path efficiency
ff = figure;
plot_scatter_mammals(mod_ind.brainSize(2).val, comm_ind(1).val,...
        col_points, x_lim, x_tick, [-0.6 0.1], -0.5:0.1:0);
xlabel('log_{10} Brain Volume (mm^3)')
ylabel('\rho (ED, SPE)')
title(sprintf('\\rho = %g; pval < 10^{-15}',...
    round(comm_ind(1).corr_log10BrVol(1), 2)))
set(gca, 'FontSize', 16, 'GridColor', [0 0 0],...
    'MinorGridColor', [0.5 0.5 0.5])
box on
exportgraphics(ff, 'ed_spe_brvol_withav.png', 'Resolution', 300)

% negative search information
ff = figure;
plot_scatter_mammals(mod_ind.brainSize(2).val, comm_ind(2).val,...
        col_points, x_lim, x_tick, [-0.64 -0.05], -0.6:0.1:-0.1);
xlabel('log_{10} Brain Volume (mm^3)')
ylabel('\rho (ED, NSI)')
indnotnan = find(not(isnan(comm_ind(2).val)));
[rho, pval] = corr(mod_ind.brainSize(2).val(indnotnan), comm_ind(2).val(indnotnan),...
    'type','Spearman');
title(sprintf('\\rho = %g; pval < 10^{-15}', round(rho ,2)))
set(gca, 'FontSize', 16, 'GridColor', [0 0 0], 'MinorGridColor', [0.5 0.5 0.5])
box on
exportgraphics(ff, 'ed_nsi_brvol_withav.png', 'Resolution', 300)

% communicability
ff = figure;
plot_scatter_mammals(mod_ind.brainSize(2).val, comm_ind(3).val,...
        col_points, x_lim, x_tick, [-0.78 -0.15], -0.7:0.1:-0.2);
xlabel('log_{10} Brain Volume (mm^3)')
ylabel('\rho (ED, CMY)')
indnotnan = find(not(isnan(distcorr_wei(:,4))));
[rho, pval] = corr(mod_ind.brainSize(2).val(indnotnan), comm_ind(3).val(indnotnan),...
    'type','Spearman');
title(sprintf('\\rho = %g; pval < 10^{-15}', round(rho ,2)))
set(gca, 'FontSize', 16, 'GridColor', [0 0 0],...
    'MinorGridColor', [0.5 0.5 0.5])
box on
exportgraphics(ff, 'ed_cmy_brvol_withav.png', 'Resolution', 300)


% get routing vs diffusion 

X = norm_sp(not(isnan(norm_co)));
Y = norm_co(not(isnan(norm_co)));

ff = figure;
sc = scatter(norm_sp(not(isnan(norm_co))), norm_co(not(isnan(norm_co))), 50,...
    col_points(not(isnan(norm_co)),:), 'filled');
sc.MarkerEdgeColor = [0 0 0];
sc.MarkerEdgeAlpha = 0.6;
axis square
grid on
grid minor
set(gca, 'XLim', [-1 0.4], 'XTick', -0.8:0.2:0.2,...
   'YLim', [-18 -4], 'YTick', -16:2:-6, 'FontSize', 15)
hold on
order_index_notnan = order_index(not(isnan(norm_co)));
for i=1:12
    pos = find(order_index_notnan==i);
    tmpcentx = mean(X(pos));
    tmpcenty = mean(Y(pos));
    tmpcolor = colors(i, :);
    hold on
    scatter(tmpcentx, tmpcenty, 90, 'MarkerEdgeColor', [0 0 0],...
        'MarkerFaceColor', [1 1 1], 'MarkerFaceAlpha', 1, 'Marker','square')
    scatter(tmpcentx, tmpcenty, 45, 'MarkerFaceColor', tmpcolor,...
        'MarkerEdgeColor', [1 1 1], 'MarkerFaceAlpha', 0.7, 'Marker','square')
end
xlabel('SPE_{scaled}')
ylabel('CMY_{scaled}')
[rho, pval] = corr(norm_sp(not(isnan(norm_co))), norm_co(not(isnan(norm_co))),...
    'type','Spearman');
title(sprintf('\\rho = %g; pval < 10^{-7}', round(rho ,2)))
set(gca, 'FontSize', 16, 'GridColor', [0 0 0],...
    'MinorGridColor', [0.5 0.5 0.5])
box on
exportgraphics(ff, 'spe_cmy_brvol_norm.png', 'Resolution', 300)


% on the 3D space also with brain volume

ff = figure;
set(gcf, 'Position', 1.0e+03*[0.0010    0.0490    1.2800    0.5993])
hold on
scatter3(mod_ind.brainSize(2).val, norm_sp, ones(201,1)*(-16),...
    'MarkerFaceColor', [0.85 0.85 0.85], 'MarkerFaceAlpha', 0.7,...
    'MarkerEdgeColor', 'none'); %, 'MarkerEdgeAlpha', 0.1);
scatter3(mod_ind.brainSize(2).val, ones(201,1)*(0.3), norm_co,...
    'MarkerFaceColor', [0.85 0.85 0.85], 'MarkerFaceAlpha', 0.7,...
    'MarkerEdgeColor', 'none'); %, 'MarkerEdgeAlpha', 0.1);
scatter3(ones(201,1)*(6.5), norm_sp, norm_co,...
    'MarkerFaceColor', [0.85 0.85 0.85], 'MarkerFaceAlpha', 0.7,...
    'MarkerEdgeColor', 'none'); %, 'MarkerEdgeAlpha', 0.1);
scatter3(mod_ind.brainSize(2).val, norm_sp, norm_co, 80, col_points, 'filled',...
    'MarkerEdgeColor', [0 0 0], 'MarkerEdgeAlpha', 0.6)
xlabel('log_{10} Brain Volume (mm^3)')
ylabel('SPE_{scaled}')
zlabel('CMY_{scaled}')
box on
grid on
set(gca, 'FontSize', 16, 'GridColor', [0 0 0],...
    'MinorGridColor', [0.5 0.5 0.5], 'View', [-50 30],...
    'XLim', [1.5 6.5], 'XTick', x_tick, 'YLim', [-0.9 0.3], 'YTick', -0.9:0.2:0.2,...
    'ZLim', [-16 -6], 'ZTick', -16:2:-6, 'YTickLabelRotation', 0)
axis square
exportgraphics(ff, 'spe_cmy_brvol_3d.png', 'Resolution', 300)



%% FIGURE 5: modularity vs communication

ff = figure;
plot_scatter_mammals(mod_ind.brainSize(2).val, corrsp,...
        col_points, x_lim, x_tick, [0.15 0.65], 0.2:0.1:0.6);
xlabel('log_{10} Brain Volume (mm^3)')
ylabel('\rho (CC, SPE)')
[rho, pval] = corr(mod_ind.brainSize(2).val, corrsp', 'type','Spearman');
title(sprintf('\\rho = %g; pval < 10^{-11}', round(rho ,2)))
set(gca, 'FontSize', 16, 'GridColor', [0 0 0],...
    'MinorGridColor', [0.5 0.5 0.5])
box on
exportgraphics(ff, 'cc_spe_brvol.png', 'Resolution', 300)

ff = figure;
plot_scatter_mammals(mod_ind.brainSize(2).val, corrco,...
        col_points, x_lim, x_tick, [0.45 0.75], 0.3:0.1:0.7);
xlabel('log_{10} Brain Volume (mm^3)')
ylabel('\rho (CC, CMY)')
[rho, pval] = corr(mod_ind.brainSize(2).val(not(isnan(corrco))), corrco(not(isnan(corrco)))', 'type','Spearman');
title(sprintf('\\rho = %g; pval < 10^{-15}', round(rho ,2)))
set(gca, 'FontSize', 16, 'GridColor', [0 0 0],...
    'MinorGridColor', [0.5 0.5 0.5])
box on
exportgraphics(ff, 'cc_cmy_brvol.png', 'Resolution', 300)

ff = figure;
plot_scatter_mammals(mod_ind.brainSize(2).val, corrsi,...
        col_points, x_lim, x_tick, [0.4 0.7], 0.45:0.1:0.65);
xlabel('log_{10} Brain Volume (mm^3)')
ylabel('\rho (CC, NSI)')
[rho, pval] = corr(mod_ind.brainSize(2).val(not(isnan(corrsi))), corrsi(not(isnan(corrsi)))', 'type','Spearman');
title(sprintf('\\rho = %g; pval < 10^{-9}', round(rho ,2)))
set(gca, 'FontSize', 16, 'GridColor', [0 0 0], 'MinorGridColor', [0.5 0.5 0.5])
box on
exportgraphics(ff, 'cc_nsi_brvol.png', 'Resolution', 300)


