function index = interpindex(xgrid,xprime)
%returns the index within xgrid where x(i) <= xprime <= x(i+1)
%if xprime is larger than the last element in xgrid, interpindex() returns
%the last index in xgrid
%if xprime is smaller than the first element in xgrid,interpindex() returns
%the first index in xgrid

%define length and allocate space
n = length(xprime);
index = zeros(n,1);

%remove final component to avoid going out of bounds
xgrid=xgrid(1:end-1);

for i =1:n
    if xprime(i)>xgrid(end)
        index(i) = length(xgrid);
    elseif xprime(i)<xgrid(1)
        index(i) = 1;
    else
        [~,index(i)] = ind2sub(size(xgrid),find(xgrid>xprime(i),1)-1);
    end
end

end
