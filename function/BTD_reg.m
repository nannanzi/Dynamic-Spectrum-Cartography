function [A,B,C] = BTD_reg(X,L0,R,C,B,p)
%BTD_REG refine A,B,C for BTD initialization
%   Detailed explanation goes here

[I,J,K]=size(X);
if nargin<6
    p=1e-2;
end
if nargin<5 | B==0
    B = randn(J,L0*R);
end
if nargin<4 | C==0
    C = randn(K,R);
end


X1 = tens2mat(X,[],1);
X2 = tens2mat(X,[],2);
X3 = tens2mat(X,[],3);
L = L0*ones(1,R);
sumL = L0*R;
A = randn(I,sumL);
Ac = Mat2Cell(A,L);
Bc = Mat2Cell(B,L);
Cc = Mat2Cell(C,ones(1,R));

iter = 100;
lambda = p;
mu = p;
for ii = 1:iter
    if ii>1*iter/2
        lambda =1e-10;
    end
    CpB = PartKron(Cc,Bc);
    A = X1'*CpB/(CpB'*CpB + lambda*eye(sumL));
    Ac = Mat2Cell(A,L);
    CpA = PartKron(Cc,Ac);
    B = X2'*CpA/(CpA'*CpA + lambda*eye(sumL));
    Bc = Mat2Cell(B,L);
    M = [];
    for kk = 1:R
        mm = reshape(Ac{kk}*Bc{kk}',[],1);
        M = [M,mm];
    end
    C = X3'*M/(M'*M  + mu*eye(R));
    Cc = Mat2Cell(C,ones(1,R));
    
    
end


end

