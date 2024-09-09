function X = sIRLSp(M,Omega,p,eta,maxiter)
% min f_p(X)  s.t. P_Omega(X) = P_Omega(M)
%   n - size of the matrix X assumed n(1) by n(2). If n is a single integer, it
% is understood that n(1) = n(2).

%   M - true matrix
%   b - data vector of the form M(Omega)

%   maxiter - maximum number of iterations
% p -  f_p(X) 0<=p<1
%  update -rules
% W = (X'X + gamma*I)^(p/2-1)
% X_new = P_Omega'(X-s*X*W) + P_Omega(M)
% gamma = gamma/eta
% s=gamma^(1-p/2)
%
if nargin < 5 || isempty(maxiter)
    maxiter=500;
end
if nargin < 4 || isempty(eta)
    eta = 1.1;
end
if nargin < 3 || isempty(p)
    p = 0;
end



m = length(Omega);
X = zeros(size(M));
[n1,n2]=size(M);
X(Omega)= M(Omega);
X0 = X;
% gamma =  1e-1*norm(M)^2;
gamma =trace(X0'*X0)/n1;
s0 =gamma^(1-p/2);
err =  zeros(maxiter,1);
for kk=1:maxiter
    W = (X'*X+gamma*eye(size(X,2)))^(p/2-1);
    Z = X-s0*X*W;
    Z(Omega) = zeros(1,length(Omega));
    X = Z + X0;
    %     if norm(Xk((Omega))-b)/norm(b)<tol
    gamma=max(gamma/eta,1e-10);
    s0 = gamma^(1-p/2);  
end


end

