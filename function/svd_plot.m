function [ output_args ] = svd_plot( X,c )
%UNTITLED3 plot the proportion of singular value
%  
% % % X =randn(100);
if nargin<2
    c = 'r-';
end
[U,S,V]=svd(X);
s = diag(S);
s = s/sum(s);
One = ones(length(s),length(s));
O = tril(One);
z = O*s;
plot(z,c,'linewidth',1);
xlabel('The first i-th index');
ylabel('Proportion');


end

