function output = LowNrankTC(M,Omega,lambda,rho,iter,Km )
%LOWNRANKTC 3-order tensor completion
%   此处显示详细说明


[I,J,K] = size(M);
if nargin < 6 || isempty(Km)
    Km =K;
end
if nargin < 5 || isempty(iter)
    iter = 100;
end


if Km>K
    Km = K;
end
M =M (:,:,K+1-Km:K);
Omega =Omega (:,:,K+1-Km:K);

X = randn(I,J,Km);
Y1 = zeros(I,J,Km);
Y2 = zeros(I,J,Km);
Y3 = zeros(I,J,Km);

% % % err = zeros(iter,1);

for ii = 1:iter
    X1  = tens2mat(X,[],1);
    X2  = tens2mat(X,[],2);
    X3  = tens2mat(X,[],3);
    Y11  = tens2mat(Y1,[],1);
    Y22  = tens2mat(Y2,[],2);
    Y33  = tens2mat(Y3,[],3);
    Z11_mat = X1 - 1/rho*Y11;
    Z11 = mySVT(Z11_mat,1/rho);
    Z22_mat = X2 - 1/rho*Y22;
    Z1 = mat2tens(Z11,[I,J,Km], [],1);
    Z22 = mySVT(Z22_mat,1/rho);
    Z33_mat = X3 - 1/rho*Y33;
    Z2 = mat2tens(Z22,[I,J,Km], [],2);
    Z33 = mySVT(Z33_mat,1/rho);
    Z3 = mat2tens(Z33,[I,J,Km], [],3);
    
    X = ( 1/(lambda+rho*3)*Omega + 1/(3*rho)*(1-Omega) ).*(Y1+Y2+Y3 + rho*Z1 +rho*Z2+ rho*Z3 )...
        +1/(lambda + rho*3)*lambda*M ;

% % %     X = Omega + 1/3*(1-Omega).*(Y1+Y2+Y3);
    Y1 = Y1 + rho*(Z1 - X);
    Y2 = Y2 + rho*(Z2 - X);
    Y3 = Y3 + rho*(Z3 - X);
% % % %     err(ii) = frob(X-Xtrue)^2/frob(Xtrue)^2;
end
output = X;
end

