function [idxx,idxy] = interpidx2D(xgrid,ygrid,xprime,yprime)
%returns the index within xgrid where x(i) <= xprime <= x(i+1) and the
%index within ygrid where y(i) <= yprime <= y(i+1)

idxx = interpindex(xgrid,xprime);
idxy = interpindex(ygrid,yprime);

end