
function  Tensor3Dvisualization(X )
%TENSOR3DVISUALIZATION �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
[I,J,K] = size(X);
XdB = 10*log10(X);
[Xgrid,Ygrid,Zgrid] = meshgrid(1:I,1:J,1:K);
h11 = slice(Xgrid,Ygrid,Zgrid,XdB,I,1,1);
hold on
h12 = slice(Xgrid,Ygrid,Zgrid,XdB,1,J,1);
h13 = slice(Xgrid,Ygrid,Zgrid,XdB,1,1,K);
set(h11,'FaceColor','interp','EdgeColor','none')
set(h12,'FaceColor','interp','EdgeColor','none')
set(h13,'FaceColor','interp','EdgeColor','none')
colormap jet

end

