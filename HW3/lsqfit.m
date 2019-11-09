function yfit = lsqfit(x,y,n)
% Returns the fit y values of order n that fits the data in x and y

%initialize arrays and variables
A = zeros(n+1,n+1);
a = zeros(n,1);
N = length(x);

%sums from normal equations
sumx = sum(x);
sumxx = sum(x.*x);
sumy = sum(y);
sumxy = sum(x.*y);

%make array to be solved using gaussian elim
A(1,1) = N;
for i = 2

%solve
a = A\b;

%evaluate fit
yfit = coeffs(1) + coeffs(2)*x;

end

