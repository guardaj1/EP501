function b = polydiv(a,N)
%
%   Divides polynomial with coefficients a by (x-N)
%

narginchk(2,3)

%setup
n = length(a);
b = zeros(1,n);

%synthetic division algorithm eq 4.26
b(1) = a(1);
for i = 2:n
    b(i) = a(i) + N*b(i-1);
end

end