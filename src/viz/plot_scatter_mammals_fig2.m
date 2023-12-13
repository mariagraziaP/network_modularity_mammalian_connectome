function sc = plot_scatter_mammals_fig2(X, Y, colpoints, colmedians, xl, xt, yl, yt)

% figure;
sc = scatter(X, Y, 50);

sc.MarkerFaceColor = colpoints;
sc.MarkerEdgeColor = [0 0 0];
% sc.MarkerEdgeAlpha = 0.4;
sc.MarkerEdgeAlpha = 0.6;

axis square
grid on
grid minor

% yl = [min(Y)-(max(Y)-min(Y))*0.1 max(Y)+(max(Y)-min(Y))*0.1];
% yl_int = (yl(2)-yl(1))/3;

set(gca, 'XLim', xl, 'XTick', xt,...
   'YLim', yl, 'YTick', yt, 'FontSize', 16, 'MinorGridColor', [0.5 0.5 0.5])

ordind = [1;1;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;...
    2;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;...
    3;3;3;3;3;3;3;3;3;3;3;3;3;3;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;5;...
    6;6;6;7;8;8;8;8;8;9;9;9;9;9;9;9;9;9;10;10;10;10;10;10;10;11;11;11;11;11;...
    11;11;11;11;11;11;11;11;11;11;11;11;11;11;11;11;11;11;11;11;11;11;11;11;...
    11;11;11;11;11;11;11;11;11;11;11;11;11;11;12;12;12;12;12;12;12;12;12;12;...
    12;12;12;12;12;12;12;12;12;12;12;12;12;12];

hold on

for i=1:12
    pos = find(ordind==i);
    
%     tmpcentx = (min(X(pos))+max(X(pos)))/2;
%     tmpcenty = (min(Y(pos))+max(Y(pos)))/2;
    tmpcentx = mean(X(pos));
    tmpcenty = mean(Y(pos));

    tmpcolor = colmedians(pos(1),:);
    
    hold on

    scatter(tmpcentx, tmpcenty, 90, 'MarkerEdgeColor', [0 0 0],...
        'MarkerFaceColor', [1 1 1], 'MarkerFaceAlpha', 1, 'Marker','square')

    scatter(tmpcentx, tmpcenty, 45, 'MarkerFaceColor', tmpcolor,...
        'MarkerEdgeColor', [1 1 1], 'MarkerFaceAlpha', 1, 'Marker','square')

end
