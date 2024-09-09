function Y= ProPorColumn(X)
%PROPORCOLUMN if two column of X are proportional, add one column with 
%   a random vector 
Xb = ColumnNormalization(X);
[m,n] = size(X);
for ii =1:n-1
    for jj = ii+1:n
        if norm(Xb(:,ii)-Xb(:,jj))<1e-12
            X(:,jj) = X(:,jj)+ 1e-4*randn(m,1);
        end
end
Y = X;
end
