function [row, col, val] = get_triu_CC(mat)

N = size(mat,1);

t = 1;
for n=1:N-1
    M = (n+1):N;
    for m=M
        row(t,1) = n;
        col(t,1) = m;
        val(t,1) = mat(n,m);
        t = t+1;
    end
end