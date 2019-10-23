function y = funcderHW2a(x)
%
%   Function for eq 3.115 in text book
%

coeffs = [5   -60   255  -450   274];

y = polyval(coeffs, x);

end
