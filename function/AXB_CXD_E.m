function X = AXB_CXD_E(A,B,C,D,E)
%AXB_CXD_E solve X from AXB + CXD = E
%  
m = size(A,1);
n = size(B,1);
[R,S,Q1,Z1]= qz(A,C,'complex');
[U,V,Q2,Z2]= qz(B',D','complex');

F = Q1*E*Q2';
Y = zeros(m,n);

for ii=n:-1:1
    if ii ==n
        Y(:,ii)=( U(n,n)'*R + V(n,n)'*S )^-1*F(:,ii);
    else
        Lefty = U(ii,ii)'*R + V(ii,ii)'*S;
% %         sk = 0;
% %         for kk = n:-1:ii+1
% %             sk = sk + U(ii,kk)'*R*Y(:,kk) + V(ii,kk)'*S*Y(:,kk);
% %         end
        sk = R*Y(:,ii+1:n)*U(ii,ii+1:n)'+S*Y(:,ii+1:n)*V(ii,ii+1:n)';
        Righty = F(:,ii) - sk;
        Y(:,ii) = Lefty\Righty;
    end
    
end
X= Z1*Y*Z2';

X= real(X);
end

