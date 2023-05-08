function [intra_conn, inter_conn, intra_conn_av, inter_conn_av, intra_conn_norm_av] =...
    get_InterIntraModulesDensity(net, comm)

N = size(net, 1);
NC = length(unique(comm));
clust_lab = unique(comm);

intra_conn = zeros(NC, 1);
inter_conn = zeros(NC, 1);

for i=1:NC

    intra_nodes = find(comm==clust_lab(i));
    intra_conn(i) = nnz(net(intra_nodes, intra_nodes));
    intra_conn_norm(i) = nnz(net(intra_nodes, intra_nodes))/length(find(comm==clust_lab(i)));

    inter_nodes = find(not(ismember(1:N, intra_nodes)));
    inter_conn(i) = nnz(net(inter_nodes, intra_nodes)) + nnz(net(intra_nodes, inter_nodes));

end

intra_conn_av = mean(intra_conn);
inter_conn_av = mean(inter_conn);

intra_conn_norm_av = mean(intra_conn_norm);