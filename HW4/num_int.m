function I = num_int(f,x)
%performs the numerical integral using the trapezoidal rule, of the
%function f over the one dimensional domain x, assuming equal grid spacing

%allocate space
I = 0;
lx = length(x);

%get grid spacing
dx = x(2)-x(1);

for ix = 1:lx-1
    I = I + .5*dx*(f(ix)+f(ix+1));
end %for

end %fnc