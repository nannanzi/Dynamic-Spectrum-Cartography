function y = SmoothVec(x,a )
%SMOOTHVEC smooth x
%   �˴���ʾ��ϸ˵��
n = length(x);
L0 = diag(-1*ones(n-1,1),1) + eye(n);
L = L0(1:end-1,:);
y = (eye(n) + a*L'*L)\x;

end

