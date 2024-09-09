function [A,B,C] = DWCPD(Y,W,Rrank,dW,A0,B0,mu )
%DWCPD Y: IxJxK; W: IxJxK;
%   Rrank: CP-rank
% dW: input window size;
% A0,B0: initilization;
if nargin<7
    mu =1e-6;
end
[I,J,K] = size(Y);
dW = min(dW,K);
Y = Y(:,:,K+1-dW:K);
W = W(:,:,K+1-dW:K);


Y1 = tens2mat(Y,[],1);
Y2 = tens2mat(Y,[],2);
Y3 = tens2mat(Y,[],3);
W1 = tens2mat(W,[],1);
W2 = tens2mat(W,[],2);
W3 = tens2mat(W,[],3);

C = randn(dW,Rrank);
A = A0;
B = B0;

%% updateC
ckBA = kr(B,A)';
for ki = 1:dW
    
    logg = W3(:,ki)==1;
    cksparse = zeros(size(ckBA));
    cksparse(:,logg) = ckBA(:,logg);
    ck = ( cksparse*cksparse' + mu*eye(Rrank) )\( cksparse*Y3(:,ki) );
    C(ki,:) = ck';
end
%% update A;
aiCB = kr(C,B)';
for ii = 1:I
    
    logg = W1(:,ii)==1;
    aisparse = zeros(size(aiCB));
    aisparse(:,logg) = aiCB(:,logg);
    ai = ( aisparse*aisparse' + mu*eye(Rrank) )\( aisparse*Y1(:,ii) );
    A(ii,:) = ai';
end
%% update B
bjCA = kr(C,A)';
for jj = 1:J
    
    logg =W2(:,jj)==1;
    bjsparse =zeros( size(bjCA) );
    bjsparse(:,logg) = bjCA(:,logg);
    bj = ( bjsparse*bjsparse' + mu*eye(Rrank) )\( bjsparse*Y2(:,jj) );
    B(jj,:) = bj';
end

end

