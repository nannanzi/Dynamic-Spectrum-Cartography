
function  Tensor3Dvisualization(X )
%TENSOR3DVISUALIZATION 此处显示有关此函数的摘要
%   此处显示详细说明
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

