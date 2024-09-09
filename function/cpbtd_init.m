function [A,B,C] = cpbtd_init(Y,Z,Q,L)
%CPBTD_INIT Y : I1 J K
%   Z: I J1 K
I = size(Z,1);
R = length(L);
U = ll1(Y,L);
Bc = {};
Cc = {};
for ii = 1:R
    Bc{ii} = U{ii}{2};
    Cc{ii} = U{ii}{3};
end
B = Cell2Mat(Bc);
C = Cell2Mat(Cc);
Z3 = tens2mat(Z,[],3);
V = Z3/(C');
Ac = {};
for ii =1:R
    vr = reshape(V(:,ii),I,[]);
    Ac{ii} = vr/(Bc{ii}'*Q');
end

A = Cell2Mat(Ac);




end

