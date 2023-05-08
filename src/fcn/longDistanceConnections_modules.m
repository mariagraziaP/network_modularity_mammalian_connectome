function [inter_mod_LD, rho, pval] = longDistanceConnections_modules(net,  coords, comm)


N = size(net, 1);

ind_conn = find(triu(net,1));
[row, col] = ind2sub([N N], ind_conn);

ED = EuclideanDistance([coords(row,1) coords(col,1)],...
    [coords(row,2) coords(col,2)],...
    [coords(row,3) coords(col,3)]);

[nod1, nod2] = ind2sub([N N], find(ED > prctile(ED, 95)));

count = 0;
for i=1:length(nod1)
    
    if comm(nod1(i)) ~= comm(nod2(i))
        count = count+1;
    end

end

inter_mod_LD = (count/length(nod1))/length(unique(comm));

