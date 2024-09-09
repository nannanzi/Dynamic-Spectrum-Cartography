%% exact CPD with factor C has size R*R
clc; clear all;
I = 8;
J = 8;
R = 10;

A = randn(I,R);
B = randn(J,R);
C = randn(R,R);

X = zeros(I,J,R);
for rr = 1:R
     arr = A(:,rr);
     brr = B(:,rr);
     AoutprodB = outprod(arr,brr);
     X = X + outprod(AoutprodB,C(:,rr));
    %      X{nn} = 
end

Qred = createQred(X);
[U,S,V] = svd(Qred);
basisR = V(:,end-R+1:end);

for rr = 1:R
    basisrr = basisR(:,rr);
    basisMat{rr} = BasisRetrieve(basisrr,R);
    validation{rr} = Qred*basisR(:,rr); 
end

%R2X validation
R2X = createR2X(X);
for rr = 1:R
    validationR2X{rr} = R2X*reshape(basisMat{rr},[],1);
end
%SD technique


basisTens = zeros(R,R,R);
for rr = 1:R
    basisTens(:,:,rr) = basisMat{rr};
end

[AA,BB,Q,Z,VV,WW] = qz(basisMat{1},basisMat{2});

Cerr = cpderr(C,VV);
%SD technique to retrieve C
