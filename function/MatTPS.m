function Mtps = MatTPS( index,M )
%MATTPS observed index
%   M size M*N with M(index) observed and M(~index)=0

[m,n] = size(M);
[Xgrid,Ygrid] = meshgrid([1:n],[1:m]);

Mtps  = TPS( Xgrid(index),Ygrid(index),M(index),Xgrid,Ygrid,1e-5);
end

