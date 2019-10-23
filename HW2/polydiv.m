function b = polydiv(a,N)
%
%   Divides polynomial with coefficients a by (x-N)
%

%setup
n = length(a);
b = zeros(1,n);

%synthetic division algorithm eq 4.26
b(n) = a(n);
for i = n-1:-1:1
    b(i) = a(i) + N*b(i+1);
end

end