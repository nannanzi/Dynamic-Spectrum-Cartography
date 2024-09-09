% create one of the basis in W, i.e., ker(R2(X)) \cap range(\piS)
function basismat = BasisRetrieve(basisvec,R)
    %create an upper diagonal first
    topR = basisvec(1:R);
    upperdiag = diag(topR);
    
    idxstart = R+1;
    for rr = 1:R
        rowlength = R - rr;
        entryidx = idxstart:(idxstart+rowlength-1);
        upperdiag(rr,rr+1:end) = basisvec(entryidx);
        idxstart = idxstart + rowlength;
    end
    
    for rr = 1:R
        rowrr = upperdiag(rr,:);
        upperdiag(:,rr) = conj(rowrr');
    end
    basismat = upperdiag;
end