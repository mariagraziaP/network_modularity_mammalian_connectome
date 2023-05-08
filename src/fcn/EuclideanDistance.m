function d = EuclideanDistance(X,Y,Z)

% X,Y,Z are matrices nPoints*2. calcola la distanza lungo la seconda
% dimensione

d = zeros(size(X, 1), 1);

for i=1:size(X,1)
    d(i) = sqrt((X(i,1)-X(i,2)).^2 + (Y(i,1)-Y(i,2)).^2 + (Z(i,1)-Z(i,2)).^2);
end