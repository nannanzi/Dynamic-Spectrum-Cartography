function Z = PartKronMat(Amat,Bmat)
%PARTKRON  Acell = {A_1,...,A_R}
%   Bcell = {B_1,...,B_R}
%  Acell o Bcell = [A_1 x B1,...., A_R x B_R]
Z = [];
col_num = size(Amat,2);
for ii = 1:col_num
    c = kron(Amat(:,ii),Bmat(:,ii));
    Z = [Z,c];
end

end
