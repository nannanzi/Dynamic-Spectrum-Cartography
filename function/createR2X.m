%create matrix R2X
%R2X matrix has columns Vec(C2(X(:,:,r1) +
%C2(X(:,:,r2)))-C2(X(:,:,r1))-C2(X(:,:,r2)));
function R2X = createR2X(X)
    [numrows,numcols,numfibs] = size(X);
    R = numfibs;
    R2X = [];
    for rr = 1:R
        matXr1 = reshape(X(:,:,rr),[numrows,numcols]);
        for rrr = 1:R
            matXr2 = reshape(X(:,:,rrr),[numrows,numcols]);
            R2Xmat = Compound(matXr1 + matXr2) - Compound(matXr1) - Compound(matXr2);
            R2Xcol = reshape(R2Xmat,[],1);
            R2X = [R2X,R2Xcol];
        end
    end
end