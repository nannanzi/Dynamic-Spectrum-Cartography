function Y = mySVT(X,r)
% SVT X = USV;
% Y = UD(S)V
%  D:
f = @(x,r) sign(x).*( ((abs(x)-r)>=0).*(abs(x)-r)   );
[U,S,V] = svd(X,'econ');
Y = U*f(S,r)*V';

end

