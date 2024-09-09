function seq = NsK(N,K )
%NSK from 1:N choose K equal-spaced number 
%   1, m+1,2*m+1,..., K*m+1, (K+1)*m+1<N
m = (N-1)/(K+1);
s = round([1:m:N]');
seq = s(2:end-1);

end

