function [finterp] = bilin_interp(x,y,f,yi,xi)

%define lengths and allocate space
n = length(xi);
m = length(yi);
finterp = zeros(m,n);

%get four points to interpolate with
[idxx,idxy] = interpidx2D(y,x,xi,yi);

%loop over all xi's and yi'x
for i = 1:n
    for j = 1:m
        %get grid point values
        x0 = x(idxx(i));
        x1 = x(idxx(i)+1);
        y0 = y(idxy(j));
        y1 = y(idxy(j)+1);
        
        if idxy(j) == 1
            idxy(j) = 2;
        end
        
        %get f values at each grid point
        f00 = f(idxy(j),idxx(i));
        f01  = f(idxy(j)-1,idxx(i));
        f10 = f(idxy(j),idxx(i)+1);
        f11 = f(idxy(j)-1,idxx(i)+1);
        
        
        %interpolate in x
        fL = f00 + ((f10-f00)/(x1-x0)) * (xi(i)-x0);
        fU = f01 + ((f11-f01)/(x1-x0)) * (xi(i)-x0);
        
        %interpolate in y
        finterp(j,i) = fL + ((fU-fL)/(y1-y0))*(yi(j)-y0);
        %finterp = flip(finterp);
    end
end