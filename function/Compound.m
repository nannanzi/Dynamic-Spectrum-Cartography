%% returning the 2nd-order compound matrix of a given matrix
function C = Compound(X)
    [nrows,ncols] = size(X);
    
    comrows = nchoosek(1:nrows,2);
    comcols = nchoosek(1:ncols,2);
    numrows = nchoosek(nrows,2);
    numcols = nchoosek(ncols,2);
    C = zeros(numrows,numcols);
    for rr = 1:numrows
        for cc = 1:numcols
            sel_row = comrows(rr,:);
            sel_col = comcols(cc,:);
            detsubmat = det(X([sel_row],[sel_col]));
            C(rr,cc) = detsubmat;
        end
    end
end