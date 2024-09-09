function Z= ColumnSumTo1( X )
%COLUMNSUMTO1 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
[m n] = size(X);
Z = zeros(m,n);
for nn=1:n
    L1 = sum(abs(X(:,nn)));
    if L1<1e-14
        Z(:,nn) = X(:,nn);
    else
        Z(:,nn) = X(:,nn)/L1;
    end
end

end

