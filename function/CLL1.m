function C = CLL1(X,L0,R )
%CLL1 X = LL1(A,B,C)
%
%%
L0=6;
% % R=3;
[I,J,K]=size(X);
Lvec = L0*ones(R,1);
sumL = L0*R;
X1 = tens2mat(X,1,[]);
X2 = tens2mat(X,2,[]);
X3 = tens2mat(X,3,[]);
for rr=1:R
    Ac{rr} = randn(I,L0);
    Bc{rr} = randn(J,L0);
    Cc{rr} = randn(K,1);
end
A = Cell2Mat(Ac);
B = Cell2Mat(Bc);
C = Cell2Mat(Cc);

iter = 50;
ferr = zeros(iter,1);
mu = 1e-8;
for ii=1:iter
    
    
    
    MA = PartKron(Cc,Bc);
    A = X1*MA/(MA'*MA + mu*eye(sumL));
%      A = ColumnNormalization(A);
    MB = PartKron(Cc,Ac);
    B = X2*MB/(MB'*MB + mu*eye(sumL));
    
    Ac = Mat2Cell(A,Lvec);
    Bc = Mat2Cell(B,Lvec);
% % %     for rr=1:R
% % %         [QB,RB] = qr(Bc{rr},0);
% % %         Bc{rr} = QB;
% % %     end
% % %     B = Cell2Mat(Bc);
    S = zeros(I*J,R);
    for rr = 1:R
        Sc{rr} = Ac{rr}*Bc{rr}';
        
        S(:,rr) = Sc{rr}(:);
    end
    C = X3*S/(S'*S+ mu*eye(R));
%     C = ColumnNormalization(C);
%     C = ColumnPositive(C);
    Cc = Mat2Cell(C,ones(R,1));
    
    
    ferr(ii) = frob(X3-C*S')^2/frob(X3)^2;
end
semilogy(ferr)
ALSerr = cpderr(Ctrue,C)

U0{1} = randn(I,sumL);
U0{2} = randn(J,sumL);
U0{3} = randn(K,R);
U = ll1_nls(X,U0,Lvec);
C = U{3};
LL1err = cpderr(Ctrue,C)
%%
end

