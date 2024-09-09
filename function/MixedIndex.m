function [Index,Rind]  =  MixedIndex( I,J,r,rho )
%DANJIINDEX 此处显示有关此函数的摘要
%   r: number of yourenji
IM = zeros(I,J);
I1 = 2*r+1;
delta = (I-1)/(I1-1);
indexI1 =round(1:delta:delta*(I1-1)+1);
index1 = indexI1(1:end-1); % m rows
Rind = index1;
IM(index1,1:J-5)=1;
temp = 0;
for rr = 1:r
    rowindex = [index1(2*rr-1):index1(2*rr)];
    cind =J-5;
    IM(rowindex,cind)=1;
    temp = mod(temp+1,2);
end

randNum = round(rho*I*J);
randIndex = randperm(I*J,randNum);
IM(randIndex) = 1;
Index = find(IM==1);

end
% % % % %%
% % % % for ii = 2:30
% % % %     [index,~] = DanjiIndex( I,J,ii );
% % % %     samplingratio = length(index)/I/J
% % % % end
