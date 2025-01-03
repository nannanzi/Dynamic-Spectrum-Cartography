function [A,B,C] = BatchCPD(Y,W,Rrank,iter,lambda,Km,A0,B0,C0)
%LOWCPRANKTC input matrix Y = Xtrue.*W
%   W 0-1 observation matrix
%  rank: CPD rank
mu =1e-6;
svalue = 1e-16;
[I,J,K] = size(Y);



if nargin < 7 || isempty(Km)
    A0 = randn(I,Rrank);
B0 = randn(J,Rrank);
end

if nargin < 6 || isempty(Km)
    Km =K;
end
if nargin < 5 || isempty(lambda)
    lambda =1;
end
if nargin < 4 || isempty(iter)
    iter=100;
end
% % if nargin < 3 || isempty(rank)
% %     iter=10;
% % end
if Km>K
    Km = K;
end
Y =Y (:,:,K+1-Km:K);
W =W (:,:,K+1-Km:K);
for kk = 1:Km
    Y(:,:,kk)=lambda^(Km-kk)*Y(:,:,kk);
end
Y1 = tens2mat(Y,[],1);
Y2 = tens2mat(Y,[],2);
Y3 = tens2mat(Y,[],3);
W1 = tens2mat(W,[],1);
W2 = tens2mat(W,[],2);
W3 = tens2mat(W,[],3);


A = A0;
B = B0;
C = C0;
% % % A = randn(I,Rrank);
% % % B = randn(J,Rrank);
% % U = cpd(Y,Rrank);
% % A = U{1};
% % B = U{2};
% % C = U{3};
C = zeros(Km,Rrank);
C(1,:) = C0(2,:);

for iitt = 1:iter
    %% updateC
    krAB = kr(B,A);
    for ki = 1:Km
        SampleIndexk = find(W3(:,ki));
        PartkrABk = krAB(SampleIndexk,:);
        yomega = Y3(SampleIndexk,ki);
%         if iitt == 1
%             ck = (PartkrABk'*PartkrABk + mu*eye(Rrank) )\(PartkrABk'*yomega);
%         else
        if ki >= 2
            ck = C(ki-1,:);
        else 
            ck = C(ki,:);
        end

        LCk = 1/norm(PartkrABk' * PartkrABk);
        gCk = (mu * ck + (ck * PartkrABk' - yomega') * PartkrABk);
        ck = ck - rand*1/LCk * gCk;
%         end
        C(ki,:) = ck;
    end
   % C(C<svalue) = svalue;
    %% update A;
    aiCB = kr(C,B)';
    for ii = 1:I
        
        logg = W1(:,ii)==1;
        aisparse = zeros(size(aiCB));
        aisparse(:,logg) = aiCB(:,logg);
        ai = ( aisparse*aisparse' + mu*eye(Rrank) )\( aisparse*Y1(:,ii) );
        A(ii,:) = ai';
    end
   % A(A<svalue) = svalue;
    %% update B
    bjCA = kr(C,A)';
    for jj = 1:J
        
        logg =W2(:,jj)==1;
        bjsparse =zeros( size(bjCA) );
        bjsparse(:,logg) = bjCA(:,logg);
        bj = ( bjsparse*bjsparse' + mu*eye(Rrank) )\( bjsparse*Y2(:,jj) );
        B(jj,:) = bj';
    end
   % B(B<svalue) = svalue;
end


end

