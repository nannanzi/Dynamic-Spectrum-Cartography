function Ind = NsK_2D(I,J,K )
%NSK_2D from IxJ choose K equal-spaced number 
%   
Mat = zeros(I,J);
Ks = ceil(sqrt(K));
row = NsK(I,Ks);
col = NsK(J,Ks);
Mat(row,col)=1;
Ind1 = find(Mat(:)==1);
rind = sort( randperm(Ks*Ks,K));
Ind = Ind1(rind);
% NewMat = zeros(I,J);
% NewMat(Ind)=1;
% pcolor(NewMat)
end

