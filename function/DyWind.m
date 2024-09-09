function z= DyWind(x,t )
%DYWIND x: vector, t: time index
%   sum( x(t-z+1:t) )>2
z = t;
zt = zeros(t,1);
if sum(x(1:t))>=2
    for ii=1:t
        zt(ii) = sum(x(t-ii+1:t));
    end
    ind = min(find(zt>=2));
    z = min(ind);
end

end

