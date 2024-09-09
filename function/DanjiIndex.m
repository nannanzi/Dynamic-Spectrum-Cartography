function [Index,Rind ] =  DanjiIndex( I,J,r )
%DANJIINDEX 此处显示有关此函数的摘要
%   此处显示详细说明
IM = zeros(I,J);
I1 = r+2;
indexI1 =1:floor((I-1)/(I1-1)):floor((I-1)/(I1-1))*I1-1;
    index1 = indexI1(2:end-1); % m rows
    Rind = index1;
IM(index1,:)=1;
temp = 0;
for rr = 1:r-1
    rowindex = [index1(rr):index1(rr+1)];
    cind = (temp==0)*J + (temp==1)*1;
    IM(rowindex,cind)=1;
    temp = mod(temp+1,2);
end


Index = find(IM==1);

end
% % % % %%
% % % % for ii = 2:30
% % % %     [index,~] = DanjiIndex( I,J,ii );
% % % %     samplingratio = length(index)/I/J
% % % % end
