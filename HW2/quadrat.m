function [x1,x2] = quadrat(coeffs)
%
%   Solves the quadratic formula for function with coefficients specified
%   in coeffs
%

x1 = (-coeffs(2) + sqrt(coeffs(2)^2 - 4*coeffs(1)*coeffs(3)))/(2*coeffs(1));
x2 = (-coeffs(2) - sqrt(coeffs(2)^2 - 4*coeffs(1)*coeffs(3)))/(2*coeffs(1));

end