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

%make A array (coeffs = A\b)
count = 0;
for ir = 1:n+1
    for ic = 1:n+1
        A(ir,ic) = sum(x.^(count));
        count = count+1;
    end
    count = ir;
end

%make b array
b = zeros(n+1,1);

for i = 1:n+1
    b(i) = sum(y.*x.^(i-1));
end


%solve systems of equations
addpath linear_algebra

[Amod,ord] = Gauss_elim(A,b);
coeffs = backsub(Amod(ord,:));

rmpath linear_algebra

%evaluate fits
yfit = zeros(length(x),1);
for j = 1:n+1
    yfit = coeffs(i)*x.^(j-1) + yfit;
end




