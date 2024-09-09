function [distidx,Permute_S] = Permutation(S2,S1)
  
    R = length(S1(1,:));
    M = length(S1(:,1));
    Permute_S = zeros(M,R);
    distidx = zeros(1,R);
    for r = 1:R
        distidxr = zeros(1,R);
        s2col = S2(:,r);
        for rr = 1:R
            s1col = S1(:,rr);
            dist21 = frob(s2col - s1col);
            distidxr(rr) = dist21;
        end
        [~,idx] = min(distidxr);
        distidx(r) = idx;
    end
    
    for rrr = 1:R
        permuteidx = distidx(rrr);
        Permute_S(:,rrr) = S2(:,permuteidx);
    end
end