function xnusqr = gof(y,yfit,sigma)
%finds the goodness-of-fit of the polynomial fit yfit when compared to the
%original data y using the uncertainty in the measurement y

%set up
nu = length(y)-1;

xnusqr = 0;
for i = length(y)
    xnusqr = xnusqr + (((y(i)-yfit(i))^2)/sigma(i)^2) ;
end
end