%create matrix R2X
%R2X matrix has columns Vec(C2(X(:,:,r1) +
%C2(X(:,:,r2)))-C2(X(:,:,r1))-C2(X(:,:,r2)));
function Qred = createQred(X)
    [numrows,numcols,numfibs] = size(X);
    R = numfibs;
    Q1 = [];
    for rr = 1:R
        matXrr = reshape(X(:,:,rr),[numrows,numcols]);
        Q1mat = Compound(matXrr + matXrr) - Compound(matXrr) - Compound(matXrr);
        Q1col = reshape(Q1mat,[],1);
        Q1 = [Q1,Q1col];
    end
    
    Q2 = [];
    comcol = nchoosek(1:R,2);
    for rrr = 1:nchoosek(R,2)
        colidx1 = comcol(rrr,1);
        colidx2 = comcol(rrr,2);
        matXr1 = reshape(X(:,:,colidx1),[numrows,numcols]);
        matXr2 = reshape(X(:,:,colidx2),[numrows,numcols]);
        Q2mat = Compound(matXr1 + matXr2) - Compound(matXr1) - Compound(matXr2);
        Q2col = reshape(Q2mat,[],1);
        Q2 = [Q2,Q2col];
    end
    
    Qred = [Q1,2*Q2];

end