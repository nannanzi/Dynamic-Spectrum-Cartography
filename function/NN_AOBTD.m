function [Ci, Si,ferr] = NN_AOBTD(X,L,R,iter,C )
%AOBTD 此处显示有关此函数的摘要
%   此处显示详细说明
[I, J, K] = size(X);
if nargin<4
    iter = 50;
end

if nargin<5
    Ci = randn(K,R).^2;
else
    Ci = C;
end
%% AOBTD
X3 = tens2mat(X,[],3);
lambda = 1e-10;

Si = randn(I*J,R);
ferr = zeros(iter+1,1);
ferr(1) = frob(X3-Si*Ci')^2;
for ii = 1:iter
% %     Si  = Si - 1/norm(Ci)^2 * (Si*(Ci'*Ci)-X3*Ci);
Si = X3*Ci*(Ci'*Ci+lambda*eye(R))^-1;

    for rr = 1:R
        [Us,Ss,Vs] = svds(reshape(Si(:,rr),I,J),L);
        Si(:,rr)= reshape(Us*Ss*Vs',[],1);
    end
    Si(Si<1e-16)=1e-16;
    % update C
    Ci = X3'*Si*(Si'*Si+lambda*eye(R))^-1;
    Ci(Ci<1e-16)=1e-16;
    ferr(ii+1) =  frob(X3-Si*Ci')^2;
end

end

