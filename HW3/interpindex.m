function index = interpindex(xgrid,xprime)
%returns the index within xgrid where x(i) <= xprime <= x(i+1)

n = length(xprime);
index = zeros(n,1);

for i =1:n
    indices = find(xgrid>xprime(i));
    index(i) = indices(1)-1;
end

end
