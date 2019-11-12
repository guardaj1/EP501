function [idxx,idxy] = interpidx2D(xgrid,ygrid,xprime,yprime)
%returns the index within xgrid where x(i) <= xprime <= x(i+1) and the
%index within ygrid where y(i) <= yprime <= y(i+1)

%find indices
for i = 1:length(xprime)
    idxx(i) = interpindex(xgrid,xprime(i));
end

for j = 1:length(yprime)
    idxy(j) = interpindex(ygrid,yprime(j));
end
end