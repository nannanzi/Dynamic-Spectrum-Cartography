function NAE = NAE_estimate(X,Xhat)
    K = length(Xhat(1,1,:));
    
    nae = zeros(1,K);
    ae = zeros(1,K);
    for k = 1:K
        
        minus = abs(X(:,:,k) - Xhat(:,:,k));
        nae(k) = sum(minus(:));
        Xhatk = abs(Xhat(:,:,k));
        ae(k) = sum(Xhatk(:));
    end
    
    NAE = sum(nae)/sum(ae);
end