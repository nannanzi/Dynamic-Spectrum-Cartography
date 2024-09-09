function z = Dist_Clocation(A,B)
%DIST_CLOCATION 此处显示有关此函数的摘要
%   此处显示详细说明
A = reshape(A,[],1);
B = reshape(B,[],1);
nA = length(A);
nB = length(B);
z = 0;
if nA==nB
    allIND = perms(1:nA);
    NN = size(allIND,1);
    Loc = zeros(NN,1);
    for ii = 1:NN 
        Loc(ii) = sqrt( sum(abs(A-B(allIND(ii,:))).^2/nA)  );
    end
    z = min(Loc);
end

end

