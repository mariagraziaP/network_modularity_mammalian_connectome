function f = plot_distributions(data, colors)



width = 0.5; % mammals
base = 1-width/2; % mammals


% figure;

for idx=1:length(data)
    
    curr_data = data{idx};
    
    nObs = size(curr_data, 1);
    
    jit = (rand(nObs,1))*width;
    jit = jit + base + idx-1;
    f{1,idx} = scatter(jit, curr_data);
    f{1,idx}.SizeData = 40;
    f{1,idx}.MarkerFaceColor = colors(idx,:,1);
%     f{1,idx}.MarkerEdgeColor = 'none';
    f{1,idx}.MarkerEdgeColor = [1 1 1];
%     f{1,idx}.MarkerEdgeAlpha = 0.4;

    if length(curr_data)>1

        hold on
    
        f{2,idx} = scatter(idx, mean(curr_data));
        f{2,idx}.SizeData = 50;
%         f{2,idx}.MarkerFaceColor = colors(idx,:,2);
%         f{2,idx}.MarkerEdgeColor = 'none';
        f{2,idx}.MarkerEdgeColor = colors(idx,:,2);
        f{2,idx}.MarkerFaceColor = [1 1 1];
        f{2,idx}.LineWidth = 1;
        f{2,idx}.MarkerFaceAlpha = 0.7;

%         if length(curr_data)>2
%             curr_std = std(curr_data);
%     
%             hold on
%     
%             f{3,idx} = line([idx idx], [mean(curr_data)-curr_std/2 mean(curr_data)+curr_std/2],...
%                 'Color', colors(idx,:,2), 'LineWidth', 1.2);
%         end
    end
    hold on
end