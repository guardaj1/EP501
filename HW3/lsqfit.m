function yfit = lsqfit(x,y)
% Returns the coefficients of the linear fit for the set of data in x and y
% of order n

%sums from normal equations
sumx = sum(x);
sumxx = sum(x.*x);
sumy = sum(y);
sumxy = sum(x.*y);

%make array to be solved using gaussian elim
A = [length(x), sumx ;...
     sumx     , sumxx];
b = [sumy,sumxy]';

%Perform gauss elim
addpath C:\Users\JulioG2793\Documents\GitHub\EP501

Amod = Gauss_elim(A,b);
coeffs = backsub(Amod);

rmpath C:\Users\JulioG2793\Documents\GitHub\EP501

%evaluate fit
yfit = coeffs(1) + coeffs(2)*x;

end

